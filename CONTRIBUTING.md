# Contributing to dyna-clust

Thank you for your interest in contributing to dyna-clust! This document provides guidelines for contributing to the project.

## How to Contribute

### Reporting Bugs

If you find a bug, please open an issue on GitHub with:

- A clear description of the bug
- Steps to reproduce the issue
- Expected behavior vs actual behavior
- Your environment (OS, R version, package versions)
- Any error messages or logs

### Suggesting Enhancements

Enhancement suggestions are welcome! Please open an issue describing:

- The motivation for the enhancement
- Detailed description of the proposed feature
- Any potential implementation approaches

### Pull Requests

1. Fork the repository
2. Create a new branch for your feature (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Test your changes thoroughly
5. Commit your changes (`git commit -m 'Add amazing feature'`)
6. Push to your branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

## Development Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/LukeLikesDirt/dyna-clust.git
   cd dyna-clust
   ```

2. Install R dependencies:
   ```r
   install.packages(c("dplyr", "readr", "stringr", "purrr", "tibble", "optparse"))
   ```

3. Install external dependencies:
   - BLAST+ (https://blast.ncbi.nlm.nih.gov/)
   - dnabarcoder (https://github.com/vuthuyduong/dnabarcoder)

## Code Style

- Follow the tidyverse style guide for R code
- Use meaningful variable and function names
- Add comments for complex logic
- Include roxygen2 documentation for all exported functions

## Testing

Before submitting a pull request:

1. Test with the example data:
   ```bash
   cd examples
   ./run_example.sh
   ```

2. Verify all functions work as expected
3. Check that documentation is up to date

## Documentation

- Update README.md if you change functionality
- Add roxygen2 comments to new functions
- Update examples if needed

## Code of Conduct

- Be respectful and inclusive
- Welcome newcomers and help them learn
- Focus on constructive feedback
- Prioritize the community and project goals

## Questions?

If you have questions about contributing, feel free to open an issue or contact the maintainers.

## License

By contributing to dyna-clust, you agree that your contributions will be licensed under the MIT License.
