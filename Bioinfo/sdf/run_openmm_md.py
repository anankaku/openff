from openff.toolkit import Molecule, ForceField
import openmm
import openmm.app as app
import openmm.unit as unit

mol = Molecule.from_file("sar.sdf")

if not mol.conformers:
    mol.generate_conformers(n_conformers=1)

ff = ForceField("openff-2.2.0.offxml")
interchange = ff.create_interchange(mol.to_topology())

system = interchange.to_openmm()
topology = interchange.topology.to_openmm()
positions = interchange.positions.to_openmm()

integrator = openmm.LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    0.002 * unit.picoseconds,
)

simulation = app.Simulation(topology, system, integrator)
simulation.context.setPositions(positions)

simulation.minimizeEnergy()

with open("topology.pdb", "w") as f:
    app.PDBFile.writeFile(topology, positions, f)

simulation.reporters.append(app.DCDReporter("traj.dcd", 1000))

simulation.reporters.append(
    app.StateDataReporter(
        "md.log",
        1000,
        step=True,
        potentialEnergy=True,
        temperature=True,
    )
)

simulation.step(500000)
