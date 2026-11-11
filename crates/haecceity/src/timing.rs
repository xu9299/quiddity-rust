//! Optional per-call diagnostics. Timers never participate in geometry or acceptance.
//! A separate observer reports active work even while the calling thread is inside a
//! long native computation. Repeated face spans are aggregated rather than streamed.

use std::cell::{Cell, RefCell};
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
    static PARALLEL: Cell<bool> = const { Cell::new(false) };
}

/// Scoped workers contribute to their caller's totals without creating observer threads.
/// Parallel face spans are siblings; their inclusive times can exceed elapsed wall time.
#[derive(Clone)]
pub(crate) struct ParallelContext(Option<Arc<Trace>>);

pub(crate) fn parallel_context() -> ParallelContext {
    ParallelContext(CURRENT.with(|current| current.borrow().clone()))
}

impl ParallelContext {
    pub(crate) fn run<R>(&self, work: impl FnOnce() -> R) -> R {
        struct Binding {
            previous: Option<Arc<Trace>>,
            parallel: bool,
        }
        impl Drop for Binding {
            fn drop(&mut self) {
                CURRENT.with(|current| current.replace(self.previous.take()));
                PARALLEL.with(|parallel| parallel.set(self.parallel));
            }
        }
        let _binding = Binding {
            previous: CURRENT.with(|current| current.replace(self.0.clone())),
            parallel: PARALLEL.with(|parallel| parallel.replace(true)),
        };
        work()
    }
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
                loop {
                    let event = match receive.recv_timeout(interval) {
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
                        Err(mpsc::RecvTimeoutError::Timeout) => {
                            let state = locked(&observed);
                            let stack = state.snapshot(Instant::now());
                            let totals = state.totals.values().cloned().collect();
                            drop(state);
                            Event {
                                event: "heartbeat",
                                elapsed_seconds: observed.started.elapsed().as_secs_f64(),
                                stack,
                                totals,
                            }
                        }
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

/// A face is visible in heartbeats; only slow completions get their own event.
pub fn face(name: &'static str, index: usize, entity: Option<u64>) -> Span {
    let operation = PARALLEL.with(|parallel| parallel.get().then_some("parallel"));
    begin(name, operation, Some(index), entity, false)
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
            .filter(|span| span.id != self.id && span.frame.operation != Some("parallel"))
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

    #[test]
    fn parallel_faces_are_siblings_and_totals_remain_in_the_callers_session() {
        let (send, receive) = mpsc::channel();
        let session = Session::with_interval(
            move |event| {
                let _ = send.send(event);
            },
            Duration::from_millis(10),
        )
        .unwrap();
        let parent = stage("mass.parallel_shell");
        let context = parallel_context();
        let barrier = std::sync::Barrier::new(3);
        thread::scope(|scope| {
            for index in 0..2 {
                let context = &context;
                let barrier = &barrier;
                scope.spawn(move || {
                    context.run(|| {
                        let face = face("mass.face_mass", index, Some(index as u64));
                        // Emit the short completion too, so the test needs no slow sleep.
                        locked(face.trace.as_ref().unwrap())
                            .stack
                            .iter_mut()
                            .find(|span| span.id == face.id)
                            .unwrap()
                            .verbose = true;
                        barrier.wait();
                        barrier.wait();
                        drop(face);
                        assert!(PARALLEL.with(Cell::get));
                    })
                });
            }
            barrier.wait();
            let heartbeat = loop {
                let event = receive.recv_timeout(Duration::from_secs(2)).unwrap();
                if event.event == "heartbeat" {
                    break event;
                }
            };
            assert_eq!(heartbeat.stack.len(), 3);
            assert!(
                heartbeat.stack[1..]
                    .iter()
                    .all(|frame| frame.operation == Some("parallel"))
            );
            barrier.wait();
        });
        assert!(!PARALLEL.with(Cell::get));
        drop(parent);
        drop(session);
        let events: Vec<_> = receive.try_iter().collect();
        let summary = events.last().unwrap();
        let completed: Vec<_> = events
            .iter()
            .filter(|event| {
                event.event == "end"
                    && event
                        .stack
                        .last()
                        .is_some_and(|frame| frame.operation == Some("parallel"))
            })
            .collect();
        assert_eq!(completed.len(), 2);
        for event in completed {
            assert_eq!(event.stack.len(), 2);
            assert_eq!(event.stack[0].name, "mass.parallel_shell");
            assert_eq!(
                event.stack[1].entity,
                event.stack[1].index.map(|index| index as u64)
            );
        }
        assert!(summary.totals.iter().any(|total| {
            total.name == "mass.face_mass"
                && total.operation == Some("parallel")
                && total.calls == 2
        }));
        assert!(CURRENT.with(|current| current.borrow().is_none()));
    }
}
