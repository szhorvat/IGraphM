This directory contains various kinds of tests for IGraph/M.

  - `RunTests.nb` runs the main test suite and some basic code quality checks. The [MicroTest](https://github.com/szhorvat/MicroTest) package is required.
  - `MakeTests.nb` contains helpers to create new unit tests.
  - `CodeInspect.nb` runs the linter on the codebase.
  - `CheckDefinitionOrdering.nb` does further code quality checks.
  - `Layouts.nb` exercises graph layouts. The output should be inspected visually.
  - `GraphDataComparisons.nb` runs time-consuming comparisons with the `GraphData` database.
  - `Editor*.wl` files contain interactive tests for `IGGraphEditor[]`.
