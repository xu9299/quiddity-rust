//! Deliver diagnostics outside the geometry thread and without altering returned JSON.
use pyo3::prelude::*;
use quiddity::kernel::timing::{Event, Session};
use serde_json::json;

fn document(event: Event) -> String {
    json!({
        "schema": "quiddity-native-timing/1",
        "event": event.event,
        "elapsed_seconds": event.elapsed_seconds,
        "stack": event.stack.into_iter().map(|frame| json!({
            "name": frame.name,
            "operation": frame.operation,
            "index": frame.index,
            "entity": frame.entity,
            "elapsed_seconds": frame.elapsed_seconds,
        })).collect::<Vec<_>>(),
        "totals": event.totals.into_iter().map(|total| json!({
            "name": total.name,
            "operation": total.operation,
            "calls": total.calls,
            "elapsed_seconds": total.elapsed_seconds,
            "max_seconds": total.max_seconds,
        })).collect::<Vec<_>>(),
    })
    .to_string()
}

pub(super) fn start(callback: Option<Py<PyAny>>) -> Option<Session> {
    let callback = callback?;
    let mut enabled = true;
    match Session::start(move |event| {
        if !enabled {
            return;
        }
        let payload = document(event);
        Python::attach(|py| {
            if let Err(error) = callback.call1(py, (payload,)) {
                // A failed observer is reported once; it cannot replace a geometry result.
                enabled = false;
                error.print(py);
            }
        });
    }) {
        Ok(session) => Some(session),
        Err(error) => {
            eprintln!("Quiddity timing observer could not start: {error}");
            None
        }
    }
}
