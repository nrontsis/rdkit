### How to apply
Apply one by one via `patch -p1 < PATCH_PATH`, in the correct order.
For any errors encountered, skip the patch while noting the files it affects.
At end copy paste any to-be-changed files from failed patches from original source.
Make sure to revert diff include "third_party" and "google3".