# Code-Style Guidelines

In addition, please pay attention to `PORTABILITY.md`.

## Basename
Basenames of files can be retrieved using `std::filesystem::path(...).filename().string()`.

## Attributes?
- `[[nodiscard]]` for a function whose return value should not be discarded
- `[[maybe_unused]]` when an attribute may be unused
