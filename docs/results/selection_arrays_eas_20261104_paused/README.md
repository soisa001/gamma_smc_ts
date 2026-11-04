# Paused selection-array continuation

Paused at the user's request. All simulation, decoder, and monitoring processes were stopped. Nothing resumes automatically.

The scheduler recorded 540 selected survivors (47 new), 1,384 additional completed proposals, no failed jobs, and 23 outstanding submitted jobs. The continuation lasted 22.93 minutes including cache replay and the reader-lock restart. Counts are provisional scheduler totals; no updated performance claims are made here.

Completed simulation receipts and original seeds are preserved. On a future authorized resume, verified completed units are reused and incomplete units are retried under a new deadline. The initial completed report remains in `../selection_arrays_eas_20261104/`.

The scheduler now retries brief Windows reader locks during atomic ledger replacement; simulation, scoring, and seed functions are unchanged. All nine focused selection-array tests passed.
