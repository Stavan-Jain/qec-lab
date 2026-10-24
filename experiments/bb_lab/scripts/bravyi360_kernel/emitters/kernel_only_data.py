"""Keep certificate definitions in the kernel while omitting unused runtime code."""
import re

def mark_kernel_data(paths):
    for path in paths:
        source = path.read_text()
        path.write_text(re.sub(r"(?m)^(\s*)(private )?def ", r"\1\2noncomputable def ", source))
