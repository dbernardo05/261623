"""Run bounded recording and restart regressions with CXX=g++-16 python3 tests/test_recording.py."""
import csv
import math
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
from collections import Counter

ROOT = Path(__file__).resolve().parents[1]
BASE = '8d5c8350fd6b602b6eea15b24088d9f67aad067e'
VARIANTS = ('uniform', 'smallworld', 'scalefree')
COMPACT = {'field_potentials.tsv', 'spike_events.tsv', 'edges.tsv', 'spikes.txt'}
STEPS, SNAP, EPOCH = 167, 154, 17
N = 6


def build(work, variant, model):
    top = work / model
    src = top / variant
    src.mkdir(parents=True)
    for name in ('main.cpp', 'layer.cpp', 'layer.h', 'neuron.cpp', 'neuron.h'):
        rel = f'{variant}/{name}'
        if model == 'corrected':
            shutil.copyfile(ROOT / rel, src / name)
        else:
            (src / name).write_bytes(subprocess.check_output(['git', 'show', f'{BASE}:{rel}'], cwd=ROOT))
    if model == 'reference':
        subprocess.run(['git', 'apply', '--whitespace=nowarn', f'--include={variant}/*',
                        ROOT / 'reference/recording.patch'], cwd=top, check=True)
    main = (src / 'main.cpp').read_text().replace('100,100', f'{N},{N}').replace('time(NULL)', '12345')
    sources = ['main.cpp', 'layer.cpp', 'neuron.cpp']
    flags = ['-DMACGREGOR_REFERENCE'] if model == 'reference' else []
    if model == 'upstream':
        main = main.replace('Param* param;', 'Param storage; Param* param = &storage;')
        p = src / 'layer.cpp'
        p.write_text(p.read_text().replace('iterations < 60000', f'iterations < {STEPS}'))
    else:
        for p in [*ROOT.glob('*.h'), ROOT / 'simulation.cpp']:
            shutil.copyfile(p, top / p.name)
        sources.append('../simulation.cpp')
    (src / 'main.cpp').write_text(main)
    exe = src / 'sim'
    subprocess.run([os.environ.get('CXX', 'g++'), '-std=c++11', '-fopenmp', '-I.', *flags, *sources,
                    '-o', exe], check=True, cwd=src)
    return exe


def rows(path):
    with path.open() as f:
        return list(csv.DictReader(f, delimiter='\t'))


def check_recordings(out, stdout, legacy, variant):
    fields = rows(out / 'field_potentials.tsv')
    events = rows(out / 'spike_events.tsv')
    spikes = [tuple(map(int, s.split())) for s in (out / 'spikes.txt').read_text().splitlines()]
    assert len(fields) == len(spikes) == STEPS
    stdout = stdout.splitlines()[-STEPS:]
    hist = Counter((int(x['step']), int(x['neuron'])) for x in events)
    assert events and max(hist.values()) == 1
    for step, (field, count) in enumerate(zip(fields, spikes)):
        sampled, tick = stdout[step].split()
        assert step == count[0] == int(field['step']) == int(tick)
        assert math.isclose(float(sampled), float(field['legacy' if legacy else 'fixed']), rel_tol=5.1e-6)
        fired = [i for i in range(N*N) if hist[step, i]]
        inh = sum(i//N % 5 == 0 and i%N % 5 == 0 for i in fired)
        assert count[1:] == (len(fired)-inh, inh)

    epochs = rows(out / 'epochs.tsv')
    cell_epochs = rows(out / 'neuron_epochs.tsv')
    assert len(cell_epochs) == len(epochs) * N*N
    release = []
    if variant == 'smallworld':
        release = [list(map(float, line.split())) for line in (out / 'releasedglutamate.txt').read_text().splitlines()]
    for idx, record in enumerate(epochs):
        start, end = idx*EPOCH, min((idx+1)*EPOCH, STEPS)
        assert (int(record['first_step']), int(record['end_step'])) == (start, end)
        assert int(record['excitatory_spikes']) == sum(s[1] for s in spikes[start:end])
        assert int(record['inhibitory_spikes']) == sum(s[2] for s in spikes[start:end])
        for name in ('legacy', 'fixed'):
            mean = sum(float(x[name]) for x in fields[start:end]) / (end-start)
            assert math.isclose(float(record['mean_'+name]), mean, abs_tol=1e-14)
        for i, summary in enumerate(cell_epochs[idx*N*N:(idx+1)*N*N]):
            assert (int(summary['first_step']), int(summary['end_step']), int(summary['neuron'])) == (start, end, i)
            assert int(summary['spikes']) == sum(hist[t, i] for t in range(start, end))
            if release and end < STEPS:
                mean = sum(release[t+1][i] for t in range(start, end)) / (end-start)
                assert math.isclose(float(summary['mean_received_glutamate']), mean, abs_tol=1e-12)

    # Decode the checkpoint independently. Recompute both electrodes on the actual post-update grid.
    data = (out / f'snapshot-{SNAP}.bin').read_bytes()
    head = struct.Struct('=8s6I18di')
    magic, endian, topology, model, height, width, step, *_ = head.unpack_from(data)
    assert (magic, endian, topology, model, height, width, step) == (
        b'MGSNAP01', 0x01020304, VARIANTS.index(variant)+1, legacy, N, N, SNAP)
    pos = head.size
    cell = struct.Struct('=7d2i')
    voltages = [cell.unpack_from(data, pos + i*cell.size)[0] for i in range(N*N)]
    pos += N*N*cell.size
    for name, decay in [('legacy', 2.0), ('fixed', 1/8)]:
        value = sum(math.exp(-decay*((i//N-N//2)**2+(i%N-N//2)**2)) / math.sqrt(8*math.pi)*e
                    for i, e in enumerate(voltages))
        assert math.isclose(value, float(fields[SNAP-1][name]), abs_tol=1e-12)
    graph = []
    for post in range(N*N):
        degree, = struct.unpack_from('=I', data, pos)
        pos += 4
        for _ in range(degree):
            r, c, weight, _, _, _, slots, used = struct.unpack_from('=2i4d2I', data, pos)
            assert used <= slots == 149
            pos += struct.calcsize('=2i4d2I') + used * struct.calcsize('=Id')
            graph.append((r*N+c, post, weight, int(bool(r%5 or c%5))))
    assert graph == [(int(e['pre_neuron']), int(e['post_neuron']), float(e['weight']), int(e['excitatory']))
                     for e in rows(out / 'edges.tsv')]


def exercise(work, variant):
    params = (ROOT / variant / 'parameters.txt').read_text().splitlines()
    params = '\n'.join('SpontRelease 0.2' if x.startswith('SpontRelease') else
                       'Neighborhood 2' if x.startswith('Neigh') else x for x in params).encode()
    def run(exe, directory, *args, success=True):
        directory.mkdir(exist_ok=True, parents=True)
        # A restart must restore the saved team even when the caller has a different default.
        threads = '2' if '--resume' in args else '1'
        result = subprocess.run([exe, *map(str, args)], cwd=directory, capture_output=True, input=params,
                                env=dict(os.environ, OMP_NUM_THREADS=threads, OMP_DYNAMIC='FALSE'), timeout=90)
        assert (result.returncode == 0) == success, (args, result.stderr)
        return result.stdout

    base = work / 'upstream' / variant / 'results'
    published = run(build(work, variant, 'upstream'), base)
    for model in ('reference', 'corrected'):
        exe = build(work, variant, model)
        out = exe.parent / 'results'
        plain = exe.parent / 'default'
        run(exe, plain, '--steps', STEPS)
        assert {p.name for p in plain.iterdir()} == COMPACT
        for name in ('weightmap.txt', 'potential.txt'):
            (plain/name).write_bytes(b'existing dense output')
        run(exe, plain, '--steps', STEPS)
        assert all((plain/name).read_bytes() == b'existing dense output' for name in ('weightmap.txt', 'potential.txt'))
        stdout = run(exe, out, '--steps', STEPS, '--per-neuron-output', '--epoch-steps', EPOCH, '--snapshot-steps', SNAP)
        check_recordings(out, stdout, model == 'reference', variant)
        if model == 'reference':
            assert published == stdout
            assert all((out/p).read_bytes() == (plain/p).read_bytes() for p in COMPACT)
            for old in base.glob('*.txt'):
                a, b = old.read_bytes(), (out/old.name).read_bytes()
                if variant == 'uniform' and old.name in ('glutamatemap.txt', 'initialglutamatemap.txt', 'weightmap.txt'):
                    # The deposited uniform source comments out these three writes.
                    assert not a and b.strip()
                    continue
                if variant == 'smallworld' and old.name == 'gaba.txt':
                    a, b = a.split(b'\n', 1)[1], b.split(b'\n', 1)[1]
                assert a == b, old.name
        saved = out / f'snapshot-{SNAP}.bin'
        prefix = f'from-{SNAP}-'
        run(exe, out, '--steps', STEPS, '--resume', saved, '--epoch-steps', EPOCH, '--per-neuron-output')
        for name in ['field_potentials.tsv', 'spike_events.tsv', 'epochs.tsv', 'neuron_epochs.tsv']:
            key, limit = ('end_step', SNAP) if 'epoch' in name else ('step', SNAP-1)
            assert rows(out/(prefix+name)) == [x for x in rows(out/name) if int(x[key]) > limit], name
        assert (out/(prefix+'edges.tsv')).read_bytes() == (out/'edges.tsv').read_bytes()
        assert (out/(prefix+'spikes.txt')).read_bytes().splitlines() == (out/'spikes.txt').read_bytes().splitlines()[SNAP:]
        if variant == 'smallworld':
            for name in ('potential', 'gaba', 'kconductance', 'releasedglutamate', 'availableglutamate', 'releasableglutamate', 'spikepatterns'):
                filename = name + '.txt'
                assert (out/(prefix+filename)).read_bytes().splitlines() == (out/filename).read_bytes().splitlines()[SNAP:]
        invalid = exe.parent / 'invalid'
        damaged = exe.parent / 'truncated.bin'
        damaged.write_bytes(saved.read_bytes()[:-1])
        run(exe, invalid, '--resume', damaged, success=False)
        run(exe, invalid, '--epoch-steps', '0', success=False)
        if model == 'corrected':
            run(exe, invalid, '--resume', work/'reference'/variant/'results'/saved.name, success=False)
        assert not list(invalid.iterdir())
        print(f'{variant} {model} passed', flush=True)


if __name__ == '__main__':
    with tempfile.TemporaryDirectory(prefix='macgregor-recording-test-') as temp:
        for variant in VARIANTS:
            exercise(Path(temp), variant)
