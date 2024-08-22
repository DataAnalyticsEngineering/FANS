# Contributing to FANS

Contributions to FANS are most welcome! Please refer to the steps below for more details.

## Changelog

We maintain a `CHANGELOG.md` where all major changes and contributions are entered.

## How to contribute

1. **Fork and Clone**: Fork the repository on GitHub and clone your fork locally.

    ```bash
    git clone https://github.com/your-username/FANS.git
    cd FANS
    ```

2. **Create a Branch**: Create a branch for your work, using a descriptive name.

    ```bash
    git checkout -b feature/my-feature
    ```

3. **Make Changes**: Implement your changes, adhering to the [Code Style Guidelines](#code-style-guidelines).

4. **Write Tests**: Ensure new features or bug fixes are covered by tests.

5. **Commit and Push**: Commit your changes with a clear message, then push to your fork.

    ```bash
    git add .
    git commit -m "Describe your changes"
    git push origin feature/my-feature
    ```

6. **Create a Pull Request**: Open a pull request to the `develop` branch. Include relevant details, such as the issue being fixed or the feature being added.

### Code Style Guidelines

- **C++ Standard**: Use C++17 or later.
- **Indentation**: 4 spaces, no tabs.
- **Naming**:
  - Functions: `camelCase`
  - Classes: `PascalCase`
  - Variables: `snake_case`
  - Constants: `ALL_CAPS`
- **Documentation**: Use Doxygen-style comments.

### Branching and Merging

- **`main`**: Latest stable release.
- **`develop`**: Active development. Base your feature branches off `develop`.
- **Feature branches**: Branch off `develop` and submit pull requests back to `develop`.
- **Release branches**: Merged into `main` for new releases.
