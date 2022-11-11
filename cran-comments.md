* This is another update to the 0.2.0 submission
- Should actually fix the clang-UBSAN compilation error
- The error was due to an Eigen bug in arrayed indexing
- The solution was to unroll the loops within the Stan model.

* This is an update from 0.1.0 to 0.2.0

- Fixes the compilation error (clang-UBSAN)
- Added predict and fitted methods for lmmelsm objects
- Fixed Documentation cross-reference mishap
- Now using inherits() instead of class() == character for class checking
