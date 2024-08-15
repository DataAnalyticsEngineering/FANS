# Docker
We provide a set of docker images for different use cases on our [Dockerhub profile](https://hub.docker.com/u/unistuttgartdae):
- **fans-ci**: Contains the minimum tools to build FANS (including dev packages of dependencies with the required headers), but does not include FANS itself. Meant for a CI workflow.
- **fans-dev**: Based upon fans-ci, but offers a non-root user (`develop`) and handling of UID and GID to not mess up permissions when volume mounting into the container. Meant as an quick to setup build environment for FANS.

Both images are built for linux/amd64 and linux/arm64 as well as for the three most recent Ubuntu LTS versions (focal, jammy, noble). The Ubuntu version can be selected through tags, e.g. `fans-dev:focal`; `noble` is equivalent to the `latest` tag. The architecture is selected automatically depending on your host platform.

## Set up a Container
Set up a development container with your current working directory (in there, use `git clone` to obtain the latest FANS version) mounted into it. You need to have [Docker Desktop](https://www.docker.com/products/docker-desktop/) installed on your machine.

First, clone FANS:
```bash
git clone https://github.tik.uni-stuttgart.de/DAE/FANS.git
cd FANS
```
Then we create the container using our `fans-dev` image.

### In a Linux, MacOS or Windows Subsystem for Linux (WSL) Shell
```bash
docker create --name fans-dev -it \
  -e HOST_UID=$(id -u) \
  -e HOST_GID=$(id -g) \
  -v /etc/localtime:/etc/localtime:ro \
  -v /etc/timezone:/etc/timezone:ro \
  -v $PWD/:/FANS/ \
  unistuttgartdae/fans-dev:latest
```
The `-e` options provide the entrypoint script of the container with your host user ID and GID, such that the user ID and GID inside the container can be adapted to match yours. This is done to not mess up file permissions in the mounted volumes. The two volume mounts of `/etc/localtime` and `/etc/timezone` are required to have the host date and time inside the container.

### In Windows PowerShell
Using PowerShell is not recommended since it only has limited support of file permissions and completely ignores file ownership in the WSL->Container direction.
```PS
docker create --name fans-dev -it `
  --env HOST_UID=1000 `
  --env HOST_GID=1000 `
  --env TZ=Europe/Berlin `
  --volume ${PWD}:/FANS/ `
  unistuttgartdae/fans-dev
```

## Working with the Container
The following workflow is suggested: You would work on the code as usual on your host; and only to build and run FANS you would attach to the container:
```bash
docker start fans-dev
docker attach fans-dev

cd /FANS
mkdir build
cd build
cmake ..
cmake --build . -j

cd ../test
./FANS
./run_tests.sh
cat nohup_test_*.log
```
For convenience we added some basic utilities to our `fans-dev` image including `htop`, `vim` and `python`.

### Attaching Visual Studio Code
You can attach VS Code to the newly created container in order to actually work inside the container. This has the benefit that IntelliSense and other static analysis tools have access to all the headers of FANS' dependencies which would not be possible when developing on the host and only using the container for building FANS.

To attach VS Code you need to install the `Remomte Development Extension Pack` and the `Docker` Extension. Then open the Docker menu, right click our newly created `fans-dev` container and select "Start" (if not running already) and then "Attach Visual Studio Code".

After attaching VS Code you unfortunately are user `root` in VS Code due to the way the UID and GID mapping is implemented: The container starts as root, executes the entrypoint script which changes UID and GID and only then drops privileges using `gosu`. VS Code though skips the entrypoint script and thus doesn't switch to the non-root user `develop`. You however can do so manually by typing `gosu develop bash` in your terminal sessions inside VS Code.

For further reading and alternative approaches like a full DevContainer setup have a look at
- [Developing inside a Container](https://code.visualstudio.com/docs/devcontainers/containers)
- [Attach to a running Container](https://code.visualstudio.com/docs/devcontainers/attach-container)
- [Specifying the default container user](https://code.visualstudio.com/remote/advancedcontainers/add-nonroot-user#_specifying-the-default-container-user)

### Calling Containerized FANS from the Host
By building inside the container, FANS is linked against the container's libs and therefore must run inside the container. After attaching to the container you can then continue to use FANS as described in the main [README](../README.md#usage). Just remember that any input and output files need to visible to the container and thus must lie somewhere inside the mounted volumes.

Special care has to be taken if you need to use FANS within scripts on the host, as Docker's interactive mode (`-i`) is not suitable in this case. Instead you need to use `docker exec`. One basically replaces the original `FANS` call by `docker exec -u develop -w /FANS/test fans-dev [original call]`. For example in conjunction with nohup:
```bash
docker start fans-dev
nohup /usr/bin/time -v docker exec -u develop -w /FANS/test fans-dev [original call] &
docker stop fans-dev
```