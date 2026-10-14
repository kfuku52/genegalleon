#!/usr/bin/env python3
"""Compare native HDF5 probability values without depending on storage filters.

Run with the qualified optional ANNEVO Python environment (numpy and h5py).
Groups, attributes, dataset shapes/dtypes and every value bit are included.
"""
import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np


def attributes(obj):
    result = {}
    for name, value in sorted(obj.attrs.items()):
        array = np.asarray(value)
        if array.dtype.hasobject:
            items = array.ravel().tolist()
            if not all(isinstance(x, (str, bytes)) for x in items):
                raise ValueError('Unsupported object-valued native HDF5 attribute')
            data = [{'text': x} if isinstance(x, str) else {'bytes': x.hex()} for x in items]
        else:
            data = {'bits_sha256': hashlib.sha256(array.tobytes(order='C')).hexdigest()}
        result[name] = {'dtype': array.dtype.str, 'shape': list(array.shape), 'data': data}
    return result


def signature(path):
    path = Path(path)
    before = path.stat()
    result = {'schema': 1, 'groups': {}, 'datasets': {}}
    with h5py.File(path, 'r') as handle:
        result['groups']['/'] = attributes(handle)
        def visit(name, obj):
            if isinstance(obj, h5py.Group):
                result['groups'][name] = attributes(obj)
            else:
                if obj.dtype.hasobject or obj.shape is None:
                    raise ValueError('Unsupported native HDF5 probability dataset')
                digest = hashlib.sha256()
                if not obj.shape:
                    digest.update(np.asarray(obj[()]).tobytes(order='C'))
                else:
                    per_row = obj.dtype.itemsize * max(1, int(np.prod(obj.shape[1:])))
                    block = max(1, 64 * 1024**2 // per_row)
                    for start in range(0, obj.shape[0], block):
                        digest.update(obj[start:start+block].tobytes(order='C'))
                result['datasets'][name] = dict(dtype=obj.dtype.str, shape=list(obj.shape),
                                               attributes=attributes(obj), bits_sha256=digest.hexdigest())
        handle.visititems(visit)
    after = path.stat()
    def identity(info):
        return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns
    if identity(before) != identity(after) or not result['datasets']:
        raise ValueError('Native probabilities changed while reading or contain no datasets')
    result['semantic_sha256'] = hashlib.sha256(json.dumps(result, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    result['stored_bytes'] = after.st_size
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    if args.input.resolve() == args.output.resolve():
        raise ValueError('Probability signature must not replace its input')
    value = signature(args.input)
    with args.output.open('x') as stream:
        json.dump(value, stream, sort_keys=True, indent=2)
        stream.write('\n')


if __name__ == '__main__':
    main()
