#!/usr/bin/env python3
import importlib.util, os, tempfile
from pathlib import Path
project=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('runtime',project/'validate_runtime_inputs.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
with tempfile.TemporaryDirectory() as td:
 root=Path(td); visible=root/'visible';visible.mkdir();marker=visible/'marker';marker.write_text('ok');image=root/'image';image.write_text('x');log=root/'args';app=root/'apptainer'
 app.write_text(f'''#!/usr/bin/env bash
set -euo pipefail
printf '%s\\n' "$@" > {log}
[[ $1 == exec ]]; shift; seen=false
while (($#)); do case "$1" in --cleanenv|--no-home) shift;; --bind) [[ $2 == {visible}:{visible} ]]; seen=true; shift 2;; *) image=$1; shift; break;; esac; done
$seen; test -s "$image"; exec "$@"
''');app.chmod(0o755)
 mod.exec_in(app,image,'sh','-c',f'test -s {marker}',bind_paths=[visible],label='binding regression')
 args=log.read_text().splitlines();assert '--bind' in args and f'{visible}:{visible}' in args
print('container binding regression test: PASS')
