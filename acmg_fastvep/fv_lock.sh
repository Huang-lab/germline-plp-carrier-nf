# shellcheck shell=bash
# Sourced by run_fastvep.sh and run_fastvep_batch.sh: one-writer locks.
#
# Two jobs writing the same chunk (a resubmit while the first is still
# running, two batches pointed at one --outdir) used to share one partial file:
# the second job's cleanup unlinked the first's partial, and the first then
# published the second's half-written output. Partial names are now per job,
# and a lock refuses a second writer outright.
#
# flock(1) where it exists (util-linux, as on Minerva): the kernel drops the
# lock when the job dies, however it dies. Elsewhere a mkdir lock, which a job
# killed with SIGKILL leaves behind; the error then says whose it was and how
# to remove it.

FV_LOCK_DIRS=""

# fv_lock_acquire <lock path> <fd number> <what is being locked, for messages>
fv_lock_acquire() {
    local lock="$1" fd="$2" what="$3" owner
    owner="pid=$$ host=$(hostname 2>/dev/null || echo unknown) LSB_JOBID=${LSB_JOBID:-none} since=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    if command -v flock >/dev/null 2>&1; then
        eval "exec ${fd}>>\"\$lock\""
        if ! flock -n "$fd"; then
            echo "ERROR: $what is being written by another running job (flock on $lock; holder: $(cat "$lock" 2>/dev/null || echo unknown))." >&2
            echo "       Not starting a second writer. Wait for it, or bkill it first." >&2
            return 1
        fi
        printf '%s\n' "$owner" > "$lock"
    else
        if ! mkdir "$lock.d" 2>/dev/null; then
            echo "ERROR: $what is locked: $lock.d exists (holder: $(cat "$lock.d/owner" 2>/dev/null || echo unknown))." >&2
            echo "       If that job is no longer running (bjobs / ps), the lock is stale: rm -r '$lock.d'" >&2
            return 1
        fi
        printf '%s\n' "$owner" > "$lock.d/owner"
        FV_LOCK_DIRS="$FV_LOCK_DIRS $lock.d"
    fi
}

fv_lock_release() {
    local d
    for d in $FV_LOCK_DIRS; do rm -rf "$d"; done
    FV_LOCK_DIRS=""
}
