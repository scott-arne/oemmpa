"""Run one command and report its resource usage, free of the caller's footprint.

Linux copies the forking process's resident high-water mark into the child and
records it as the child's own at ``exec``, so a build spawned straight from the
benchmark harness reports the harness's footprint -- the interpreter plus
``oemmpa``, ``duckdb`` and everything they pull in -- rather than its own. The
mark is a high-water mark, so the harness cannot shed it by freeing memory
first. What it can do is fork from somewhere lean, which is this module: an
isolated interpreter that imports nothing but the standard library, spawns the
real command, and writes the figures back as JSON.

Run it as a script, never import it::

    python -I -S benchmarks/_child_rusage.py REPORT COMMAND [ARG ...]

``REPORT`` is a path this writes a JSON object to: ``status`` (the command's
exit status, a negative number for a signal, exactly as
``subprocess.Popen.returncode`` reports it), ``seconds`` (its wall time) and
``maxrss`` (its peak resident size in the platform's ``ru_maxrss`` unit). The
command's own streams are inherited untouched, which is why the figures travel
through a file rather than stdout. A non-zero exit from this module means the
measurement itself failed; a command that fails is a successful measurement
with a non-zero ``status``.
"""

import json
import os
import subprocess
import sys
from time import perf_counter


def main(argv):
    """Spawn ``argv[1:]``, reap it, and write its usage to ``argv[0]``.

    :param argv: Report path followed by the command's argument vector.
    :returns: Process exit status for this module: zero once the figures are
        written, whatever the command itself did.
    :raises IndexError: When no command follows the report path.
    """
    report_path, command = argv[0], argv[1:]
    start = perf_counter()
    process = subprocess.Popen(command)
    _, status, usage = os.wait4(process.pid, 0)
    seconds = perf_counter() - start
    # Marks the child reaped, so Popen neither warns that it is still running
    # nor polls its pid again, and gives the signal-aware status for free.
    process.returncode = os.waitstatus_to_exitcode(status)
    with open(report_path, "w", encoding="utf-8") as report:
        json.dump(
            {"status": process.returncode, "seconds": seconds, "maxrss": usage.ru_maxrss},
            report,
        )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
