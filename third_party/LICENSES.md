# Third-party components

PEG bundles unmodified copies of the header-only libraries below, so that it builds without installing anything. They are **not** part of PEG's own source code and remain under their own licenses. PEG itself is licensed under the GPLv3 (see `../LICENSE`).

| Component | Version | Location | License | License file |
|---|---|---|---|---|
| [Eigen](https://eigen.tuxfamily.org/) | 3.4.0 | `eigen/Eigen/` | Mozilla Public License 2.0 (MPL-2.0) | `COPYING.MPL2` |
| [Boost.Odeint](https://www.boost.org/) | 1.90.0 | `boost/boost/` | Boost Software License 1.0 (BSL-1.0) | `LICENSE_1_0.txt` |

## Eigen

Copyright (C) the Eigen authors. Eigen is licensed under the MPL-2.0 (https://www.mozilla.org/MPL/2.0/). The MPL is a file-level license: if you modify one of the bundled Eigen files, that file stays under the MPL-2.0 and its source must be made available. Individual files may carry additional license notices in their headers (e.g. files under `eigen/Eigen/src/misc/`); these are not used by PEG. Keep every notice in the headers intact.

## Boost

Copyright (C) the respective Boost authors. Boost is licensed under the Boost Software License 1.0 (https://www.boost.org/LICENSE_1_0.txt). Only the headers needed by `boost/numeric/odeint.hpp` and its dependencies are included, not the complete Boost distribution. The license text must accompany source distributions; it does not have to be reproduced for compiled binaries.

## License compatibility

Both licenses are compatible with the GPLv3 under which PEG is distributed (MPL-2.0 via its Secondary License clause, BSL-1.0 as a permissive license). The combined work, PEG together with these headers, is distributed under the GPLv3.

## Updating

The files are copied unmodified. To update, replace the content of `eigen/Eigen/` and `boost/boost/` with the new release, keep the license files, adjust the versions in the table above, and reinstall peg via `make`.
