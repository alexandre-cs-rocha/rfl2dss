# rlf2dss

Converts electrical network models from **Roseau Load Flow (RLF)** JSON format into a set of **OpenDSS** script files (`.dss`). The generated files can be compiled and solved directly in OpenDSS via the [AltDSS](https://github.com/dss-extensions/AltDSS-Python) Python interface.

---

## Overview

[Roseau Load Flow](https://roseau-load-flow.roseautechnologies.com/) is a Python library for modelling and solving distribution network power flows. OpenDSS is the industry-standard open-source distribution system simulator developed by EPRI.

`rlf2dss` bridges these two tools: it reads a solved RLF network and writes the equivalent OpenDSS circuit as an organised collection of `.dss` files that follow common Brazilian distribution utility naming conventions (MT = medium voltage / média tensão, BT = low voltage / baixa tensão).

---

## Requirements

- Python ≥ 3.10
- A valid **Roseau Load Flow licence key** stored in `licence_key.txt` (same directory as the script)

Install Python dependencies:

```bash
pip install -r requirements.txt
```

---

## File Structure

After running the converter, the output directory will contain one `.dss` file per component category:

| File | Contents |
|---|---|
| `Master_<name>.dss` | Entry point — clears, redirects all files, solves |
| `CircuitoMT_<name>.dss` | Voltage source (infinite bus) definition |
| `CodCondutor_<name>.dss` | Conductor codes (LineCode elements) |
| `SegmentosBT_<name>.dss` | LV line segment definitions |
| `SegmentosMT_<name>.dss` | MV line segment definitions *(placeholder)* |
| `ChavesBT_<name>.dss` | LV switch definitions |
| `ChavesMT_<name>.dss` | MV switch definitions *(placeholder)* |
| `TransformadorMTMTMTBT_<name>.dss` | Transformer definitions |
| `CargasBT_<name>.dss` | LV load definitions |
| `CargasMT_<name>.dss` | MV load definitions *(placeholder)* |
| `CargasBTNT_<name>.dss` | LV non-technical loss loads *(placeholder)* |
| `CargasMTNT_<name>.dss` | MV non-technical loss loads *(placeholder)* |
| `CurvacargaBT_<name>.dss` | LV load shapes (LoadShape) *(placeholder)* |
| `CurvacargaMT_<name>.dss` | MV load shapes (LoadShape) *(placeholder)* |
| `ReguladorMT_<name>.dss` | MV voltage regulators *(placeholder)* |
| `RamaisBT_<name>.dss` | LV service drops *(placeholder)* |
| `Medidores_<name>.dss` | Energy meters *(placeholder)* |
| `Tensoesbase_<name>.dss` | Base voltage list (`Set voltagebases`) |

---

## Usage

### 1. Basic conversion

```python
from pathlib import Path
from rlf2dss import rlf2dss
import roseau.load_flow as rlf

# Activate licence
rlf.activate_license(Path('licence_key.txt').read_text())

json_file  = Path('networks/MyFeeder.json')
output_dir = Path('output/MyFeeder')

converter = rlf2dss(json_file, 'MyFeeder')
converter.generate_dss_files(output_dir)
```

### 2. Compile and solve with AltDSS

```python
import altdss

master = output_dir / f'Master_{converter.name}.dss'
altdss.altdss.Commands(f'Compile ({master})')
altdss.altdss.Solution.Solve()
```

### 3. Run the bundled example

```bash
python rlf2dss.py
```

This assumes `11_MVFeeder0725.json` and `licence_key.txt` are in the same directory as the script.

---

## Technical Notes

### Kron Reduction

When a line's primitive impedance matrix is 4×4 (three phases + neutral), the neutral conductor is eliminated via **Kron reduction** before computing the equivalent sequence parameters:

```
Z_abc = Z_xyz − (Z_xn · Z_nx) / Z_nn
```

### Positive-Sequence Impedance

The converter currently derives a single positive-sequence impedance for each line code from the phase matrix average self (`z_s`) and mutual (`z_m`) terms:

```
z1 = z_s − z_m
```

For studies requiring fully unbalanced line models, the `Rmatrix` / `Xmatrix` form should be used instead.

### Load Modelling

Loads are mapped to OpenDSS model codes as follows:

| RLF type | OpenDSS model | Description |
|---|---|---|
| `power` | 1 | Constant power (P+jQ) |
| `impedance` | 2 | Constant impedance |
| `current` | 5 | Constant current |

Balanced three-phase loads are exported as a single multi-phase element. Unbalanced loads are split into individual single-phase elements.

---

## Extending Placeholder Files

Several files are generated as empty stubs (marked *placeholder* in the table above). To add content, implement the logic inside the corresponding private method in the `rlf2dss` class:

- `_generate_cargas_mt` — MV loads
- `_generate_curva_carga_bt` / `_generate_curva_carga_mt` — LoadShape time series
- `_generate_regulador_mt` — voltage regulators
- `_generate_ramais_bt` — LV service drops
- `_generate_medidores` — energy meters

---

## Licence

This project requires a valid **Roseau Load Flow** licence. See [Roseau Technologies](https://www.roseautechnologies.com/) for details.
