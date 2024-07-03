#!/bin/bash --login

# Abort script at first error, when a command exits with non-zero status (except in until or while loops, if-tests, list constructs)
set -e

### workaround to fix permissions in mounted volumes ###
# This is necessary because the user in the container has a different UID and GID than the user on the host.
# USAGE: docker run -e HOST_UID=$(id -u) -e HOST_GID=$(id -g) ...
# open issue on this topic: https://github.com/docker/roadmap/issues/398
hostgroup="hostgroup"
container_user="develop"

if [ "$(id -u -n)" = "root" ]; then
    if [ -n "$HOST_UID" ] && [ -n "$HOST_GID" ]; then
        echo "Setting UID and GID to match provided host UID and GID..."
        # echo "'id' before changes: $(id $container_user)"

        if ! getent group $hostgroup >/dev/null; then
            groupadd -o -g $HOST_GID $hostgroup
        fi

        old_group=$(id -g -n $container_user)

        if ! id -nG $container_user | grep -qw $hostgroup; then
            usermod -g $hostgroup $container_user
        fi

        if ! id -nG $container_user | grep -qw $old_group; then
            usermod -a -G $old_group $container_user
        fi

        if [ "$(id -u $container_user)" != "$HOST_UID" ]; then
            usermod -u $HOST_UID $container_user
        fi

        # echo "'id' after changes: $(id $container_user)"
    else
        echo "WARNING: Please provide HOST_UID and HOST_GID as environment variables (docker run -e)! UID and GID will not be changed. This will probably lead to permission issues with mounted volumes."
    fi
else
    echo "WARNING: Can't change UID and GID to given host UID and GID. entrypoint.sh must run as root! UID and GID will not be changed. This will probably lead to permission issues with mounted volumes."
fi

# drop privileges and execute given commands as the user $container_user
exec gosu $container_user "$@"
