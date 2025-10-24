# Guide to release new version of FANS

The developer who is releasing a new version of FANS is expected to follow this workflow:

The release of the `FANS` repository is made directly from a release branch called `FANS-v1.2.3`. This branch is primarily intended to assist other developers with testing.

1. Create a branch called `FANS-v1.2.3` from the latest commit of the `develop` branch.

2. Bump the version in the `CHANGELOG.md`, the base `CMakeLists.txt`, `pixi.toml`, and in the file `FANS_Dashboard/pyproject.toml` on the branch `FANS-v1.2.3`.

3. Assuming you have Pixi installed, run the command `pixi lock` in the repository root to update the `pixi.lock` file. Then commit and push the `FANS-v1.2.3` branch to remote.

4. [Open a Pull Request `main` <-- `FANS-v1.2.3`](https://github.com/DataAnalyticsEngineering/FANS/compare/main...main) named after the version (i.e., `Release v1.2.3`) and briefly describe the new features of the release in the PR description.

5. Once the CI runs successfully (all green ticks) and one approving review is made, merge the release PR (from `FANS-v1.2.3`) into `main` by a merge commit (**not** *squash and merge* or *rebase and merge*).

6. [Draft a new release](https://github.com/DataAnalyticsEngineering/FANS/releases/new) in the `Releases` section of the repository page in a web browser. The release tag needs to be the exact version number (i.e., `v1.2.3` or `v1.2.3rc1`, compare to [existing tags](https://github.com/DataAnalyticsEngineering/FANS/tags)). Use `@target:main`. Release title is also the version number (i.e., `v1.2.3` or `v1.2.3rc1`, compare to [existing releases](https://github.com/DataAnalyticsEngineering/FANS/tags)). Use the `Auto-generate release notes` feature.

7. Merge `main` into `develop` for synchronization of `develop`.

8. If everything is in order up to this point, then the new version can be released by hitting the "Publish release" button in your Release Draft. This will create the corresponding tag.
