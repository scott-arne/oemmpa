# Vendored cmake-openeye modules

These CMake modules are vendored from
[scott-arne/cmake-openeye](https://github.com/scott-arne/cmake-openeye) so the
project configures without a network fetch at configure time (clean machines,
CI, and corporate-proxy environments).

- **Upstream:** https://github.com/scott-arne/cmake-openeye
- **Version:** v1.2.1 (commit `b78af2f55c56231160ddf1d00abaa6f41e9e81b0`)

Only the consumed modules are vendored (`FindOpenEye.cmake`,
`FindOpenEyePython.cmake`, `OpenEyeSWIG.cmake`) plus the upstream `LICENSE`.

## Refreshing

To update to a newer upstream release, replace the module files and `LICENSE`
from the corresponding tag and update the version reference above:

```bash
git clone --depth 1 --branch <tag> https://github.com/scott-arne/cmake-openeye /tmp/cmake-openeye
cp /tmp/cmake-openeye/{FindOpenEye,FindOpenEyePython,OpenEyeSWIG}.cmake \
   /tmp/cmake-openeye/LICENSE \
   cmake/cmake-openeye/
```
