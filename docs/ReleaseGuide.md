# Guide to release new version of FANS

The developer who is releasing a new version of FANS is expected to follow this work flow:

The release of the `FANS` repository is made directly from a release branch called `FANS-v1.2.3`. This branch is mainly needed to help other developers with testing.

1. Create a branch called `FANS-v1.2.3` from the latest commit of the `develop` branch.

2. If it is a real release, [open a Pull Request `main` <-- `FANS-v1.2.3`](https://github.com/DataAnalyticsEngineering/FANS/compare/main...main) named after the version (i.e. `Release v1.2.3`) and briefly describe the new features of the release in the PR description.

3. Bump the version in the `CHANGELOG.md` and the base `CMakeLists.txt` on the branch `FANS-v1.2.3`.

4. [Draft a new release](https://github.com/DataAnalyticsEngineering/FANS/releases/new) in the `Releases` section of the repository page in a web browser. The release tag needs to be the exact version number (i.e.`v1.2.3` or `v1.2.3rc1`, compare to [existing tags](https://github.com/DataAnalyticsEngineering/FANS/tags)). Use `@target:main`. Release title is also the version number (i.e. `v1.2.3` or `v1.2.3rc1`, compare to [existing releases](https://github.com/DataAnalyticsEngineering/FANS/tags)).

    * *Note:* If it is a pre-release then the option *This is a pre-release* needs to be selected at the bottom of the page. Use `@target:FANS-v1.2.3` for a pre-release, since we will never merge a pre-release into `main`.
    * Use the `Auto-generate release notes` feature.

    a) If a pre-release is made: Directly hit the "Publish release" button in your Release Draft.

    b) If this is a "real" release: As soon as one approving review is made, merge the release PR (from `FANS-v1.2.3`) into `main`.

5. Merge `main` into `develop` for synchronization of `develop`.

6. If everything is in order up to this point then the new version can be released by hitting the "Publish release" button in your Release Draft. This will create the corresponding tag.
