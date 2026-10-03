"""Release only files whose paths AND bytes match the reviewed inventory.

This utility cannot approve changed or added files. It never imports the analysis,
opens empirical inputs, walks another folder for data, or adds generated outputs.
"""
import argparse
import hashlib
import json
from pathlib import Path
import stat
import zipfile

ROOT=Path(__file__).resolve().parents[1]

def check(root=ROOT):
    manifest=root/'release_inventory.json'
    if not manifest.is_file():raise ValueError('No reviewed release inventory.')
    inventory=json.loads(manifest.read_text(encoding='utf-8'))
    expected=inventory['files']
    if 'release_inventory.json' in expected:raise ValueError('Manifest must not hash itself.')
    actual=set()
    for p in root.rglob('*'):
        if p.is_symlink() or (hasattr(p,'is_junction') and p.is_junction()):raise ValueError('Links and junctions cannot enter a release.')
        if p.is_file() and '__pycache__' not in p.parts:
            actual.add(p.relative_to(root).as_posix())
    if actual != set(expected)|{'release_inventory.json'}:
        raise ValueError('Package has missing or unreviewed files: '+str(sorted(actual ^ (set(expected)|{'release_inventory.json'}))))
    for relative,record in expected.items():
        p=root/relative
        if p.resolve().parent!=root.resolve() and root.resolve() not in p.resolve().parents:raise ValueError('Path escapes package.')
        blob=p.read_bytes()
        if hashlib.sha256(blob).hexdigest()!=record['sha256'] or len(blob)!=record['bytes']:
            raise ValueError('Reviewed bytes changed: '+relative)
    return inventory

def build(output):
    output=Path(output).resolve()
    if output.exists():raise ValueError('Will not overwrite an existing archive.')
    if ROOT==output or ROOT in output.parents:raise ValueError('Archive must be outside the package directory.')
    inventory=check();output.parent.mkdir(parents=True,exist_ok=True)
    names=sorted([*inventory['files'],'release_inventory.json'])
    with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for relative in names:
            blob=(ROOT/relative).read_bytes()
            if relative!='release_inventory.json' and hashlib.sha256(blob).hexdigest()!=inventory['files'][relative]['sha256']:
                raise ValueError('File changed during archive construction.')
            info=zipfile.ZipInfo('paper-c-replication/'+relative,date_time=(2026,10,3,0,0,0))
            info.compress_type=zipfile.ZIP_DEFLATED
            info.external_attr=(stat.S_IFREG|0o644)<<16
            archive.writestr(info,blob)
    with zipfile.ZipFile(output) as archive:
        if archive.testzip() is not None:raise ValueError('Archive CRC check failed.')
        if archive.namelist()!=['paper-c-replication/'+x for x in names]:raise ValueError('Archive inventory mismatch.')
        for relative in inventory['files']:
            if hashlib.sha256(archive.read('paper-c-replication/'+relative)).hexdigest()!=inventory['files'][relative]['sha256']:
                raise ValueError('Archive bytes mismatch.')
    print('Verified archive:',output)
    print('SHA256:',hashlib.sha256(output.read_bytes()).hexdigest())

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    if args.output:build(args.output)
    else:check();print('Reviewed inventory and bytes match.')
