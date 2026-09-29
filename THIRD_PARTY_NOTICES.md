# Third-party notices

## PoseLib polynomial helpers

Upstream: https://github.com/PoseLib/PoseLib

Frozen revision: `1758bd6d20f607c276388d28b5ba88f3651e85e9`.

The following two files are copied without modification from the locally frozen
PoseLib dependency:

- `third_party/poselib/PoseLib/misc/univariate.h`
- `third_party/poselib/PoseLib/misc/univariate.cc`

Copyright (c) 2020, Viktor Larsson. BSD-3-Clause. The full license is preserved in
`third_party/poselib/LICENSE`, and both source files retain their copyright headers.
Their checksums are recorded in `SOURCE_SNAPSHOT.json`. They provide exceptional
polynomial root fallbacks; the PoseLib P3P solver implementations are not bundled.

## Eigen

Upstream: https://eigen.tuxfamily.org/

Eigen is required externally and is not bundled in this repository. Refer to the
license notices in the selected upstream Eigen distribution. This project does
not replace or relicense those notices.
