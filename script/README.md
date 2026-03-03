# redc DAG Tool

`redc_dag.py` builds dependency graphs for redc scheduling analysis.

## Template DAG

Build a semantic DAG for modeled `redc_core`:

```bash
python script/redc_dag.py template \
  --streams 2 \
  --schedule round-robin \
  --output-json out/redc_template.json \
  --output-dot out/redc_template.dot
```

## Asm DAG

Parse GNU `objdump` text into a def-use DAG:

```bash
objdump -d -C --start-address=0x427750 --stop-address=0x429140 \
  out/build/linux/benchmark/blst_benchmarks > /tmp/redc_k2.s

python script/redc_dag.py asm \
  --input /tmp/redc_k2.s \
  --output-json out/redc_asm.json \
  --output-dot out/redc_asm.dot
```

Optional address range filter:

```bash
python script/redc_dag.py asm \
  --input /tmp/redc_k2.s \
  --start-addr 0x427750 \
  --end-addr 0x429140
```
