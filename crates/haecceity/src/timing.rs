//! Optional per-call diagnostics. Timers never participate in geometry or acceptance.
//! A separate observer reports active work even while the calling thread is inside a
//! long native computation. Repeated face spans are aggregated rather than streamed.

use std::cell::RefCell;
use std::collections::BTreeMap;
use std::marker::PhantomData;
use std::rc::Rc;
use std::sync::{Arc, Mutex, MutexGuard, mpsc};
use std::thread::{self, JoinHandle};
use std::time::{Duration, Instant};

#[derive(Clone, Debug)]
pub struct Frame {
    pub name: &'static str,
    pub operation: Option<&'static str>,
    pub index: Option<usize>,
    pub entity: Option<u64>,
    pub elapsed_seconds: f64,
}

#[derive(Clone, Debug)]
pub struct Total {
    pub name: &'static str,
    pub operation: Option<&'static str>,
    pub calls: u64,
    /// Inclusive elapsed time: nested totals must not be added together.
    pub elapsed_seconds: f64,
    pub max_seconds: f64,
}

#[derive(Clone, Debug)]
pub struct Event {
    pub event: &'static str,
    pub elapsed_seconds: f64,
    pub stack: Vec<Frame>,
    pub totals: Vec<Total>,
}

struct Active {
    id: u64,
    frame: Frame,
    started: Instant,
    verbose: bool,
}

#[derive(Default)]
struct State {
    next_id: u64,
    stack: Vec<Active>,
    totals: BTreeMap<(&'static str, Option<&'static str>), Total>,
}

impl State {
    fn snapshot(&self, now: Instant) -> Vec<Frame> {
        self.stack
            .iter()
            .map(|span| Frame {
                elapsed_seconds: now.duration_since(span.started).as_secs_f64(),
                ..span.frame.clone()
            })
            .collect()
    }
}

enum Message {
    Event(Event),
    Stop,
}

struct Trace {
    started: Instant,
    state: Mutex<State>,
    send: mpsc::Sender<Message>,
}

fn locked(trace: &Trace) -> MutexGuard<'_, State> {
    trace
        .state
        .lock()
        .unwrap_or_else(std::sync::PoisonError::into_inner)
}

thread_local! {
    static CURRENT: RefCell<Option<Arc<Trace>>> = const { RefCell::new(None) };
}

/// The caller owns the session and drops it before leaving the detached native call.
/// Separate calling threads have separate sessions. No process-global hook is installed.
pub struct Session {
    trace: Arc<Trace>,
    previous: Option<Arc<Trace>>,
    worker: Option<JoinHandle<()>>,
    _thread_bound: PhantomData<Rc<()>>,
}

impl Session {
    pub fn start(sink: impl FnMut(Event) + Send + 'static) -> std::io::Result<Self> {
        Self::with_interval(sink, Duration::from_secs(5))
    }

    fn with_interval(
        mut sink: impl FnMut(Event) + Send + 'static,
        interval: Duration,
    ) -> std::io::Result<Self> {
        let (send, receive) = mpsc::channel();
        let trace = Arc::new(Trace {
            started: Instant::now(),
            state: Mutex::new(State::default()),
            send,
        });
        let observed = Arc::clone(&trace);
        let worker = thread::Builder::new()
            .name("quiddity-timing".into())
            .spawn(move || {
                let mut next_heartbeat = Instant::now() + interval;
                loop {
                    let now = Instant::now();
                    if now >= next_heartbeat {
                        let state = locked(&observed);
                        let stack = state.snapshot(Instant::now());
                        let totals = state.totals.values().cloned().collect();
                        drop(state);
                        sink(Event {
                            event: "heartbeat",
                            elapsed_seconds: observed.started.elapsed().as_secs_f64(),
                            stack,
                            totals,
                        });
                        next_heartbeat = Instant::now() + interval;
                    }
                    // Completed spans must not postpone the next periodic observation.
                    // Check the deadline before receiving too, even with a nonempty queue.
                    let wait = next_heartbeat.saturating_duration_since(Instant::now());
                    let event = match receive.recv_timeout(wait) {
                        Ok(Message::Event(event)) => event,
                        Ok(Message::Stop) | Err(mpsc::RecvTimeoutError::Disconnected) => {
                            let totals = locked(&observed).totals.values().cloned().collect();
                            sink(Event {
                                event: "summary",
                                elapsed_seconds: observed.started.elapsed().as_secs_f64(),
                                stack: Vec::new(),
                                totals,
                            });
                            break;
                        }
                        Err(mpsc::RecvTimeoutError::Timeout) => continue,
                    };
                    sink(event);
                }
            })?;
        let previous = CURRENT.with(|current| current.replace(Some(Arc::clone(&trace))));
        Ok(Self {
            trace,
            previous,
            worker: Some(worker),
            _thread_bound: PhantomData,
        })
    }
}

impl Drop for Session {
    fn drop(&mut self) {
        CURRENT.with(|current| current.replace(self.previous.take()));
        let _ = self.trace.send.send(Message::Stop);
        if let Some(worker) = self.worker.take() {
            let _ = worker.join();
        }
    }
}

pub struct Span {
    trace: Option<Arc<Trace>>,
    id: u64,
    _thread_bound: PhantomData<Rc<()>>,
}

fn begin(
    name: &'static str,
    operation: Option<&'static str>,
    index: Option<usize>,
    entity: Option<u64>,
    verbose: bool,
) -> Span {
    let trace = CURRENT.with(|current| current.borrow().clone());
    let mut id = 0;
    if let Some(trace) = &trace {
        let now = Instant::now();
        let mut state = locked(trace);
        id = state.next_id;
        state.next_id += 1;
        state.stack.push(Active {
            id,
            frame: Frame {
                name,
                operation,
                index,
                entity,
                elapsed_seconds: 0.0,
            },
            started: now,
            verbose,
        });
        if verbose {
            let event = Event {
                event: "start",
                elapsed_seconds: trace.started.elapsed().as_secs_f64(),
                stack: state.snapshot(now),
                totals: Vec::new(),
            };
            let _ = trace.send.send(Message::Event(event));
        }
    }
    Span {
        trace,
        id,
        _thread_bound: PhantomData,
    }
}

/// A coarse stage emits both boundaries, including when its computation returns an error.
pub fn stage(name: &'static str) -> Span {
    begin(name, None, None, None, true)
}

pub fn family(name: &'static str, operation: &'static str) -> Span {
    begin(name, Some(operation), None, None, true)
}

/// A repeated operation contributes totals without emitting every boundary.
pub fn operation(name: &'static str) -> Span {
    begin(name, None, None, None, false)
}

/// A face is visible in heartbeats; only slow completions get their own event.
pub fn face(name: &'static str, index: usize, entity: Option<u64>) -> Span {
    begin(name, None, Some(index), entity, false)
}

impl Drop for Span {
    fn drop(&mut self) {
        let Some(trace) = &self.trace else {
            return;
        };
        let now = Instant::now();
        let mut state = locked(trace);
        let Some(position) = state.stack.iter().rposition(|span| span.id == self.id) else {
            return;
        };
        // A completion reports its own face last, never a concurrently active sibling.
        let mut stack: Vec<_> = state
            .stack
            .iter()
            .filter(|span| span.id != self.id)
            .map(|span| Frame {
                elapsed_seconds: now.duration_since(span.started).as_secs_f64(),
                ..span.frame.clone()
            })
            .collect();
        let span = state.stack.remove(position);
        let elapsed = now.duration_since(span.started).as_secs_f64();
        stack.push(Frame {
            elapsed_seconds: elapsed,
            ..span.frame.clone()
        });
        let total = state
            .totals
            .entry((span.frame.name, span.frame.operation))
            .or_insert(Total {
                name: span.frame.name,
                operation: span.frame.operation,
                calls: 0,
                elapsed_seconds: 0.0,
                max_seconds: 0.0,
            });
        total.calls += 1;
        total.elapsed_seconds += elapsed;
        total.max_seconds = total.max_seconds.max(elapsed);
        if span.verbose || elapsed >= 0.5 {
            let event = Event {
                event: if thread::panicking() { "unwind" } else { "end" },
                elapsed_seconds: trace.started.elapsed().as_secs_f64(),
                stack,
                totals: Vec::new(),
            };
            let _ = trace.send.send(Message::Event(event));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn queued_boundaries_do_not_starve_periodic_heartbeats() {
        let (send, receive) = mpsc::channel();
        let session = Session::with_interval(
            move |event| {
                let _ = send.send(event);
                // A busy sink ensures the pending boundaries outlast the interval.
                thread::sleep(Duration::from_millis(2));
            },
            Duration::from_millis(10),
        )
        .unwrap();
        let ongoing = operation("ongoing");
        for _ in 0..30 {
            drop(stage("short"));
        }
        drop(ongoing);
        drop(session);
        let events: Vec<_> = receive.try_iter().collect();
        let heartbeat = events
            .iter()
            .position(|event| event.event == "heartbeat")
            .expect("boundaries queued before Stop must not suppress heartbeats");
        assert!(
            events[heartbeat + 1..]
                .iter()
                .any(|event| event.event == "end")
        );
        let summary = events.last().unwrap();
        assert_eq!(summary.event, "summary");
        assert!(
            summary
                .totals
                .iter()
                .any(|total| total.name == "short" && total.calls == 30)
        );
    }

    #[test]
    fn heartbeat_observes_an_unfinished_face_and_totals_are_drained() {
        let (send, receive) = mpsc::channel();
        let session = Session::with_interval(
            move |event| {
                let _ = send.send(event);
            },
            Duration::from_millis(10),
        )
        .unwrap();
        let stage = stage("step.read");
        let face = face("mass.face_mass", 9, Some(1234));
        let heartbeat = loop {
            let event = receive.recv_timeout(Duration::from_secs(2)).unwrap();
            if event.event == "heartbeat" {
                break event;
            }
        };
        assert_eq!(heartbeat.stack.len(), 2);
        assert_eq!(heartbeat.stack[1].entity, Some(1234));
        assert_eq!(heartbeat.stack[1].index, Some(9));
        drop(face);
        drop(stage);
        drop(session);
        let events: Vec<_> = receive.try_iter().collect();
        let summary = events.last().unwrap();
        assert_eq!(summary.event, "summary");
        assert!(
            summary
                .totals
                .iter()
                .any(|t| t.name == "mass.face_mass" && t.calls == 1)
        );
        assert!(CURRENT.with(|current| current.borrow().is_none()));
    }

    #[test]
    fn simultaneous_sessions_do_not_mix_face_identities() {
        let handles: Vec<_> = (0..2)
            .map(|index| {
                thread::spawn(move || {
                    let (send, receive) = mpsc::channel();
                    let session = Session::with_interval(
                        move |event| {
                            let _ = send.send(event);
                        },
                        Duration::from_millis(10),
                    )
                    .unwrap();
                    let face = face("mass.face_moments", index, Some(100 + index as u64));
                    let event = receive.recv_timeout(Duration::from_secs(2)).unwrap();
                    assert_eq!(event.stack[0].index, Some(index));
                    drop(face);
                    drop(session);
                })
            })
            .collect();
        for handle in handles {
            handle.join().unwrap();
        }
    }
}
