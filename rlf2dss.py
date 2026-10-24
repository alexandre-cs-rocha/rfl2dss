"""
rlf2dss.py
==========
Converts electrical network models from Roseau Load Flow (RLF) JSON format
into OpenDSS script files (.dss). The generated scripts can be compiled and
solved directly by OpenDSS via the AltDSS Python interface.

Typical workflow:
    1. Load a network from a Roseau Load Flow JSON file.
    2. Instantiate `rlf2dss` with the file path and a circuit name.
    3. Call `generate_dss_files()` to produce all .dss scripts in an output directory.
    4. Optionally compile and solve the circuit with AltDSS.
    
Created by: Alexandre Rocha
"""

import json
from pathlib import Path
from enum import Enum

import roseau.load_flow as rlf
import numpy as np
import altdss

# ---------------------------------------------------------------------------
# Symmetrical component constants
# ---------------------------------------------------------------------------

# Alpha: unit phasor rotated 120° counter-clockwise (used in sequence analysis)
alpha = 1 * np.cos(120 * np.pi / 180) + 1j * np.sin(120 * np.pi / 180)

# Fortescue transformation matrix [A]: converts sequence (0,+,-) to phase (a,b,c)
#   [V_abc] = [A] * [V_012]
A = np.array([[1, 1, 1], [1, alpha ** 2, alpha], [1, alpha, alpha ** 2]])


# ---------------------------------------------------------------------------
# Helper function
# ---------------------------------------------------------------------------

def reducao_kron(z_linha):
    """
    Apply Kron's reduction to eliminate the neutral conductor from a 4x4
    primitive impedance matrix, returning a 3x3 phase impedance matrix.

    Kron reduction formula:
        Z_abc = Z_xyz - (Z_xn * Z_nx) / Z_nn

    Parameters
    ----------
    z_linha : np.ndarray, shape (4, 4)
        Full primitive impedance matrix including phase conductors (a, b, c)
        and neutral (n).

    Returns
    -------
    z_abc : np.ndarray, shape (3, 3)
        Reduced impedance matrix with the neutral conductor eliminated.
    """
    z_xyz = z_linha[:3, :3]          # Phase-to-phase sub-matrix (3x3)
    z_nn = z_linha[3, 3]              # Neutral self-impedance (scalar)

    # Outer product of the phase-neutral coupling columns/rows
    z_abc = z_xyz - 1 / z_nn * np.dot(z_linha[:3, 3], z_linha[3, :3])

    return z_abc


# ---------------------------------------------------------------------------
# Enumerations
# ---------------------------------------------------------------------------

class Phases(Enum):
    """Maps RLF phase labels to OpenDSS bus node numbers."""
    a = 1
    b = 2
    c = 3
    n = 4


class Load_Model(Enum):
    """
    Maps RLF load types to OpenDSS load model codes.
    OpenDSS model codes:
        1 = Constant Power (P+jQ)
        2 = Constant Impedance (Z)
        5 = Constant Current (I)
    """
    power = 1
    impedance = 2
    current = 5


class Units(Enum):
    kilometer = 'km'
    meters = 'm'
    miles = 'mi'
    feet = 'Ft'
    

# ---------------------------------------------------------------------------
# Main converter class
# ---------------------------------------------------------------------------

class rlf2dss:
    """
    Converts a Roseau Load Flow electrical network to a set of OpenDSS (.dss)
    script files.

    The network is read from a RLF JSON file, solved once to obtain operating
    conditions, and then exported as a collection of organised .dss files that
    mirror common Brazilian distribution utility conventions (MT = medium
    voltage, BT = low voltage, NT = non-technical losses, etc.).

    Parameters
    ----------
    file : Path
        Path to the Roseau Load Flow JSON network file.
    name : str
        Label used as a suffix on every generated .dss filename and as the
        OpenDSS circuit name.
    pos_seq : bool
        Define if the linecod is going to be modeled as positive impedance
        or as Rmatrix and Xmatrix
    """

    def __init__(self, file: Path, name: str, pos_seq: bool = False) -> None:
        self.file = file
        self.name = name
        self.pos_seq = pos_seq

        # Header comment inserted at the top of every generated .dss file
        self.base_text = '!Criação do Arquivo\n'

        # Load and solve the network so that operating-point data is available
        # (e.g. nominal voltages, ampacities) when building the DSS files.
        self.en = rlf.ElectricalNetwork.from_json(file, include_results=False)
        self.en.solve_load_flow()

        # Ordered list of logical DSS file categories.  The Master file
        # references all other files via Redirect statements in this order.
        self.list_files = [
            'Master', 'CircuitoMT', 'CodCondutor', 'CurvacargaMT', 'CurvacargaBT',
            'ChavesMT', 'SegmentosMT', 'ReguladorMT', 'TransformadorMTMTMTBT',
            'ChavesBT', 'SegmentosBT', 'RamaisBT', 'Medidores', 'CargasMT',
            'CargasBT', 'CargasMTNT', 'CargasBTNT', 'Tensoesbase',
        ]

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def generate_dss_files(self, output_dir: Path) -> str:
        """
        Generate the complete set of OpenDSS script files in *output_dir*.

        All required subdirectories are created automatically.  One .dss file
        is produced per category listed in ``self.list_files``.

        Parameters
        ----------
        output_dir : Path
            Destination directory for the generated files.
        """
        self.output_dir = output_dir

        self._generate_cargas_bt()
        self._generate_cargas_btnt()
        self._generate_cargas_mt()
        self._generate_cargas_mtnt()
        self._generate_chaves_bt()
        self._generate_chaves_mt()
        self._generate_circuito_mt()
        self._generate_cod_condutor()
        self._generate_curva_carga_bt()
        self._generate_curva_carga_mt()
        self._generate_master()
        self._generate_medidores()
        self._generate_ramais_bt()
        self._generate_regulador_mt()
        self._generate_segmentos_bt()
        self._generate_segmentos_mt()
        self._generate_tensoesbase()
        self._generate_trafos()

    def native_json_str(self) -> str:
        """
        Serialise the loaded network back to a JSON string using RLF's own
        format (results are excluded).

        Returns
        -------
        str
            JSON representation of the network.
        """
        json_dict = self.en.to_dict(include_results=False)
        json_str = json.dumps(json_dict)
        return json_str

    def generate_json(self, json_str: str, output_dir: Path) -> None:
        """
        Parse a JSON string and write it to *output_dir/output.json* with
        pretty-printing (indent=4).

        Parameters
        ----------
        json_str : str
            Valid JSON string to write.
        output_dir : Path
            Target directory (created if it does not exist).
        """
        json_obj = json.loads(json_str)
        output_dir.mkdir(parents=True, exist_ok=True)
        with open(output_dir / 'output.json', 'w') as fp:
            json.dump(json_obj, fp, indent=4)

    # ------------------------------------------------------------------
    # Private file-generation helpers
    # ------------------------------------------------------------------

    def _generate_cargas_bt(self) -> None:
        """
        Generate the low-voltage (BT) load definitions file.

        Iterates over all loads in the network.  Balanced three-phase loads
        (where all phase powers are equal) are written as a single multi-phase
        Load element.  Unbalanced loads are split into individual single-phase
        Load elements, one per energised phase.

        Loads with zero total power are skipped.
        Connection type is inferred from the presence of a neutral conductor
        ('n' in the phase list to Wye; otherwise to Delta).
        """
        file_name = f'CargasBT_{self.name}'

        list_loads = [self.base_text]
        for load in self.en.loads.values():
            # Map RLF load type to OpenDSS model integer
            model = Load_Model[f'{load.type}'].value

            # Nominal voltage in kV (OpenDSS uses kV)
            kv = load.bus.nominal_voltage.magnitude / 1000

            # Determine OpenDSS connection type
            conn = 'Delta'
            if 'n' in load.phases:
                conn = 'Wye'
            
            # Extract complex power magnitudes for each phase
            match model:
                case 1:
                    powers = [p for p in load.powers.magnitude / 1000]

                case 2:
                    impedance = np.array([i for i in load.impedances.magnitude], dtype=np.complex128())
                    kw = kv**2 * impedance.real / np.abs(impedance)**2
                    kvar = kv**2 * impedance.imag / np.abs(impedance)**2
                    if conn == 'Delta':
                        kw *= 3
                        kvar *= 3
                    powers = kw + 1j*kvar

            if len(powers) == 3 and powers[0] == powers[1] and powers[0] == powers[2]:
                # --- Balanced load: single multi-phase element ---
                if powers[0] == 0:
                    continue
                phases = [str(Phases[f'{p}'].value) for p in load.phases]
                phases = '.'.join(phases)
                kw = powers[0].real
                pf = 1
                if np.abs(powers[0]) != 0:
                    pf = kw / np.abs(powers[0])   # Power factor = P / |S|
                    
                list_loads.append(
                    f'New "Load.{load.id}" bus1="{load.bus.id}.{phases}" '
                    f'phases={len(load.voltage_phases)} conn={conn} model={model} '
                    f'kv={kv} kw={kw} pf={pf} status=variable vmaxpu=1.5 vminpu=0.92\n'
                )

            else:
                # --- Unbalanced load: one element per non-zero phase ---
                for i, power in enumerate(powers):
                    if power == 0:
                        continue
                    phases = [str(Phases[f'{p}'].value) for p in load.voltage_phases[i]]
                    phases = '.'.join(phases)
                    kw = power.real
                    pf = kw / np.abs(power)
                    list_loads.append(
                        f'New "Load.{load.id}_{load.voltage_phases[i]}" '
                        f'bus1="{load.bus.id}.{phases}" phases=1 conn={conn} '
                        f'model={model} kv={kv} kw={kw} pf={pf} '
                        f'status=variable vmaxpu=1.5 vminpu=0.92\n'
                    )

        self._file_witer(list_loads, file_name)

    def _generate_cargas_btnt(self) -> None:
        """
        Generate the low-voltage non-technical-loss (BT-NT) load file.

        Currently outputs an empty placeholder file.  Populate with NT load
        elements when non-technical loss modelling is required.
        """
        file_name = f'CargasBTNT_{self.name}'
        list_loads = [self.base_text]
        self._file_witer(list_loads, file_name)

    def _generate_cargas_mt(self) -> None:
        """
        Generate the medium-voltage (MT) load definitions file.

        Currently outputs an empty placeholder file.
        """
        file_name = f'CargasMT_{self.name}'
        list_loads = [self.base_text]
        self._file_witer(list_loads, file_name)

    def _generate_cargas_mtnt(self) -> None:
        """
        Generate the medium-voltage non-technical-loss (MT-NT) load file.

        Currently outputs an empty placeholder file.
        """
        file_name = f'CargasMTNT_{self.name}'
        list_loads = [self.base_text]
        self._file_witer(list_loads, file_name)

    def _generate_chaves_bt(self) -> None:
        """
        Generate the low-voltage (BT) switch definitions file.

        Maps every RLF Switch element to an OpenDSS Line element with
        ``switch=y``, connecting the same buses and phases.
        """
        file_name = f'ChavesBT_{self.name}'

        list_switches = [self.base_text]
        for switch in self.en.switches.values():
            # Build phase node strings for both terminals, e.g. "1.2.3"
            phases1 = [str(Phases[f'{p}'].value) for p in switch.phases1]
            phases1 = '.'.join(phases1)
            phases2 = [str(Phases[f'{p}'].value) for p in switch.phases2]
            phases2 = '.'.join(phases2)
            list_switches.append(
                f'New "Line.{switch.id}" phases={len(switch.phases)} '
                f'bus1="{switch.bus1.id}.{phases1}" bus2="{switch.bus2.id}.{phases2}" '
                f'r1=0.001 r0=0.001 x1=0 x0=0 c1=0 c0=0 length=0.001 switch=T\n'
            )

        self._file_witer(list_switches, file_name)

    def _generate_chaves_mt(self) -> None:
        """
        Generate the medium-voltage (MT) switch definitions file.

        Currently outputs an empty placeholder file.
        """
        file_name = f'ChavesMT_{self.name}'
        list_loads = [self.base_text]
        self._file_witer(list_loads, file_name)

    def _generate_circuito_mt(self) -> None:
        """
        Generate the OpenDSS circuit (voltage source) definition file.

        Creates one ``Circuit`` element per RLF voltage source.  The source
        impedance is set to nearly zero (r1=0, x1=0.0001) to represent an
        ideal infinite bus, and the per-unit voltage is fixed at 1.0.
        """
        file_name = f'CircuitoMT_{self.name}'

        list_circuito = [self.base_text]
        for source in self.en.sources.values():
            bus = source.bus
            basekv = bus.nominal_voltage.magnitude / 1000   # Convert V to kV
            list_circuito.append(
                f'New "Circuit.{self.name}" basekv={basekv} pu=1 '
                f'bus1="{bus.id}" r1=0 x1=0.0001\n'
            )

        self._file_witer(list_circuito, file_name)

    def _generate_cod_condutor(self) -> None:
        """
        Generate the conductor code (LineCode) definitions file.

        For each line in the network, the impedance representation written to the
        file depends on the ``self.pos_seq`` flag:

        **Positive-sequence mode** (``pos_seq=True``):
        - If the impedance matrix is 4x4 (phases + neutral), Kron reduction is
        applied first to obtain a 3x3 phase matrix.
        - The positive-sequence impedance is computed from the reduced matrix:
            z_s = average of diagonal terms (self-impedance)
            z_m = average of off-diagonal terms (mutual impedance)
            z_1 = z_s - z_m
        - Written as a scalar ``r1`` / ``x1`` pair (Ω/km).
        - Suitable for balanced three-phase studies.

        **Full-matrix mode** (``pos_seq=False``):
        - The raw phase impedance matrix is used directly (no reduction applied).
        - Resistance and reactance are separated into ``Rmatrix`` and ``Xmatrix``.
        - Each matrix row is formatted as space-separated values; rows are joined
        with `` | `` to match the OpenDSS lower-triangular syntax.
        - Suitable for unbalanced analysis.
        """
        file_name = f'CodCondutor_{self.name}'
        list_codcon = [self.base_text]
        cods = set()

        for line in self.en.lines.values():
            if line.parameters.id in cods:
                continue
            cods.add(line.parameters.id)
            z_abc = line.z_line.magnitude / line.length.magnitude   # Raw impedance matrix (3x3 or 4x4), Ω/km

            if self.pos_seq == True:
                # --- Positive-sequence (scalar) representation ---

                # Apply Kron reduction if a neutral conductor is present (4x4 matrix)
                if z_abc.shape == (4, 4):
                    z_abc = reducao_kron(z_abc)

                # Average self and mutual-impedance terms across the three phases
                z_s = 1/3 * sum([z_abc[i][i] for i in range(3)])           # Self impedance
                z_m = 1/3 * sum([z_abc[0][1], z_abc[1][2], z_abc[2][0]])   # Mutual impedance
                z_1 = z_s - z_m                                             # Positive-sequence impedance

                list_codcon.append(
                    f'New "Linecode.{line.parameters.id}" nphases={len(line.phases)} '
                    f'basefreq=60 r1={z_1.real} x1={z_1.imag} units=km '
                    f'normamps={line.ampacities[0].magnitude}\n'
                )

            else:
                # --- Full-matrix (unbalanced) representation ---

                # Build the Rmatrix string: each row as space-separated floats,
                # rows separated by ' | ' (OpenDSS lower-triangular row delimiter)
                rmatrix = z_abc.real
                str_rmatrix = []
                for rline in rmatrix:
                    str_line = [str(r) for r in rline]
                    str_rmatrix.append(' '.join(str_line))
                str_r = ' | '.join(str_rmatrix)

                # Build the Xmatrix string in the same format
                xmatrix = z_abc.imag
                str_xmatrix = []
                for xline in xmatrix:
                    str_line = [str(x) for x in xline]
                    str_xmatrix.append(' '.join(str_line))
                str_x = ' | '.join(str_xmatrix)

                # Write the Linecode using the '~' continuation syntax for matrices
                list_codcon.append(
                    f'New "Linecode.{line.parameters.id}" nphases={len(line.phases)} '
                    f'basefreq=60 units=km normamps={line.ampacities[0].magnitude}\n'
                    f'~ rmatrix = [{str_r}]\n'
                    f'~ xmatrix = [{str_x}]\n'
                )

        self._file_witer(list_codcon, file_name)

    def _generate_curva_carga_bt(self) -> None:
        """
        Generate the low-voltage (BT) load curve (LoadShape) file.

        Currently outputs an empty placeholder file.
        """
        file_name = f'CurvacargaBT_{self.name}'
        list_loads = [self.base_text]
        self._file_witer(list_loads, file_name)

    def _generate_curva_carga_mt(self) -> None:
        """
        Generate the medium-voltage (MT) load curve (LoadShape) file.

        Currently outputs an empty placeholder file.
        """
        file_name = f'CurvacargaMT_{self.name}'
        list_loads = [self.base_text]
        self._file_witer(list_loads, file_name)

    def _generate_master(self) -> None:
        """
        Generate the OpenDSS master file.

        The master file:
        1. Clears any existing circuit with ``Clear``.
        2. Redirects to every other .dss file (in the canonical order defined
           by ``self.list_files``) so they are processed in the correct sequence.
        3. Sets simulation mode to *daily*, configures convergence tolerances,
           and issues a ``Solve`` command.
        """
        file_name = f'Master_{self.name}'
        file_list = ['Clear\n']

        # Add a Redirect line for each component file (skip the Master itself)
        for f in self.list_files:
            if f != 'Master':
                file_list.append(f"Redirect '{f}_{self.name}.dss'\n")

        # Simulation settings
        file_list.append('Set mode = daily\n')
        file_list.append('Set tolerance = 0.0001\n')
        file_list.append('Set maxcontroliter = 10\n')
        file_list.append('!Set algorithm = newton\n')   # Newton-Raphson (optional)
        file_list.append('!Solve mode = direct\n')      # Direct solve (optional)
        file_list.append('Solve\n')
        file_list.append('!Export meters\n')            # Meter export (optional)

        self._file_witer(file_list, file_name)

    def _generate_medidores(self) -> None:
        """
        Generate the energy meter (EnergyMeter) definitions file.

        Currently outputs an empty placeholder file.
        """
        file_name = f'Medidores_{self.name}'
        list_loads = [self.base_text]
        self._file_witer(list_loads, file_name)

    def _generate_ramais_bt(self) -> None:
        """
        Generate the low-voltage service drop (ramal BT) definitions file.

        Currently outputs an empty placeholder file.
        """
        file_name = f'RamaisBT_{self.name}'
        list_loads = [self.base_text]
        self._file_witer(list_loads, file_name)

    def _generate_regulador_mt(self) -> None:
        """
        Generate the medium-voltage (MT) voltage regulator definitions file.

        Currently outputs an empty placeholder file.
        """
        file_name = f'ReguladorMT_{self.name}'
        list_loads = [self.base_text]
        self._file_witer(list_loads, file_name)

    def _generate_segmentos_bt(self) -> None:
        """
        Generate the low-voltage (BT) line segment definitions file.

        Each RLF Line element is written as an OpenDSS ``Line`` element
        referencing the corresponding ``Linecode`` (from _generate_cod_condutor).
        Phase node strings are built for both terminals.
        """
        file_name = f'SegmentosBT_{self.name}'

        list_lines = []
        for line in self.en.lines.values():
            # Build phase-node strings, e.g. phases "a,b,c" to "1.2.3"
            phases1 = [str(Phases[f'{p}'].value) for p in line.phases1]
            phases1 = '.'.join(phases1)
            phases2 = [str(Phases[f'{p}'].value) for p in line.phases2]
            phases2 = '.'.join(phases2)

            # Map RLF length unit to OpenDSS string
            units = Units[f'{line.length.units}'].value

            list_lines.append(
                f'New "Line.{line.id}" phases={len(line.phases)} '
                f'bus1="{line.bus1.id}.{phases1}" bus2="{line.bus2.id}.{phases2}" '
                f'linecode="{line.parameters.id}" length={line.length.magnitude} '
                f'units={units}\n'
            )

        self._file_witer(list_lines, file_name)

    def _generate_segmentos_mt(self) -> None:
        """
        Generate the medium-voltage (MT) line segment definitions file.

        Currently outputs an empty placeholder file.
        """
        file_name = f'SegmentosMT_{self.name}'
        list_loads = [self.base_text]
        self._file_witer(list_loads, file_name)

    def _generate_tensoesbase(self) -> None:
        """
        Generate the base voltage settings file.

        Collects the set of unique nominal voltages (in kV) from all buses in
        the network and writes a single OpenDSS ``Set voltagebases`` command.
        OpenDSS uses this list to correctly transform the solution in per unit values.
        """
        file_name = f'Tensoesbase_{self.name}'

        tensoesbase = set()
        list_tensoes = [self.base_text]
        for bus in self.en.buses.values():
            tensoesbase.add(bus.nominal_voltage.magnitude / 1000)   # V to kV

        list_tensoes.append(f'Set voltagebases={list(tensoesbase)}\n')
        list_tensoes.append('Calcvoltagebases\n')
        self._file_witer(list_tensoes, file_name)

    def _generate_trafos(self) -> None:
        """
        Generate the transformer definitions file (MT-MT and MT-BT).

        For each RLF Transformer:
        - Winding connections are extracted from the vector group string and
          mapped to OpenDSS 'Wye' / 'Delta' labels.
        - Rated apparent power (kVA), percentage load loss (%loadloss) and
          percentage no-load loss (%noloadloss) are computed from the RLF
          parameters.
        - The tap position is taken directly from the RLF model.
        """
        file_name = f'TransformadorMTMTMTBT_{self.name}'

        list_trafos = [self.base_text]
        for trafo in self.en.transformers.values():
            # Extract individual winding connection letters from the vector group
            # (e.g. 'Dyn11' to ['D', 'y']) and map to OpenDSS labels
            wg = trafo.parameters.extract_windings(trafo.parameters.vg)[:-1]
            wg = ['Delta' if w.lower() == 'd' else 'Wye' for w in wg]

            # Build bus-node connection strings for HV and LV terminals
            conn1 = [str(Phases[f'{p}'].value) for p in trafo.bus_hv.phases]
            conn1 = '.'.join(conn1)
            conn2 = [str(Phases[f'{p}'].value) for p in trafo.bus_lv.phases]
            conn2 = '.'.join(conn2)

            # Rated apparent power in VA
            sn = trafo.sn.magnitude

            # Losses as percentages of rated power
            loadloss   = trafo.parameters.psc.magnitude / sn * 100   # Short-circuit (load) loss
            noloadloss = trafo.parameters.p0.magnitude  / sn * 100   # No-load (core) loss

            list_trafos.append(
                f'New "Transformer.{trafo.id}" phases={len(trafo.phases1)} '
                f'windings=2 '
                f'buses=["{trafo.bus_hv.id}.{conn1}" "{trafo.bus_lv.id}.{conn2}"] '
                f'conns=[{wg[0]} {wg[1]}] '
                f'kvs=[{trafo.bus_hv.nominal_voltage.magnitude / 1000} '
                f'{trafo.bus_lv.nominal_voltage.magnitude / 1000}] '
                f'tap={trafo.tap} kva={sn / 1000} '
                f'%loadloss={loadloss} %noloadloss={noloadloss}\n'
            )

        self._file_witer(list_trafos, file_name)

    # ------------------------------------------------------------------
    # Internal utility
    # ------------------------------------------------------------------

    def _file_witer(self, list_str: list, file_name: str) -> None:
        """
        Write a list of strings to a .dss file in the configured output directory.

        The output directory is created automatically if it does not exist.

        Parameters
        ----------
        list_str : list of str
            Lines to write to the file (each should include a trailing newline).
        file_name : str
            Base name of the output file (without the .dss extension).
        """
        self.output_dir.mkdir(parents=True, exist_ok=True)
        with open(self.output_dir / f'{file_name}.dss', 'w') as fp:
            fp.writelines(list_str)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

#TODO: Comparar resultados dos fluxos de carga.

if __name__ == '__main__':
    # Activate the Roseau Load Flow licence before any network operations
    rlf.activate_license(Path('license_key.txt').read_text())

    # Resolve paths relative to this script's location
    dir_path = Path(__file__).parent
    name_feeder = 'MV_Feeder_VPA01'

    # Input: RLF JSON network file
    json_file = dir_path / 'input' / 'EDPES' / 'rlf' / f'{name_feeder}.json'

    # Output: directory where all .dss files will be written
    output_dir = dir_path / 'output' / f'{name_feeder}'

    # Convert the network to OpenDSS format
    dss_writer = rlf2dss(json_file, name_feeder)
    dss_writer.generate_dss_files(output_dir)

    # Compile and solve the generated master file with AltDSS
    altdss.altdss.Commands(f'Compile ({output_dir / f"Master_{dss_writer.name}"})')
    altdss.altdss.Solution.Solve()
