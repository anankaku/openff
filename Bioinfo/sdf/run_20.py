from openff.toolkit import Molecule, ForceField
import openmm
import openmm.app as app
import openmm.unit as unit

# -----------------------
# Settings
# -----------------------
sdf_file = "sar.sdf"      # 改成你的 SDF 檔名
n_conformers = 20
steps_per_conformer = 100000   # 100000 steps x 2 fs = 0.2 ns
report_interval = 1000

# -----------------------
# Build OpenFF molecule
# -----------------------
mol = Molecule.from_file(sdf_file)

# Generate multiple starting conformers
mol.generate_conformers(n_conformers=n_conformers)

ff = ForceField("openff-2.2.0.offxml")
interchange = ff.create_interchange(mol.to_topology())

system = interchange.to_openmm()
topology = interchange.topology.to_openmm()

# -----------------------
# OpenMM setup
# -----------------------
integrator = openmm.LangevinMiddleIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    0.002 * unit.picoseconds,
)

platform = openmm.Platform.getPlatformByName("CUDA")
properties = {
    "Precision": "mixed",
    "DeviceIndex": "0",
}

simulation = app.Simulation(
    topology,
    system,
    integrator,
    platform,
    properties,
)

# Save topology once
with open("topology.pdb", "w") as f:
    app.PDBFile.writeFile(
        topology,
        mol.conformers[0].to_openmm(),
        f,
    )

# -----------------------
# Run short MD from each conformer
# -----------------------
for i, conf in enumerate(mol.conformers):
    print(f"Running conformer {i}")

    positions = conf.to_openmm()
    simulation.context.setPositions(positions)
    simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

    simulation.minimizeEnergy()

    # Clear old reporters
    simulation.reporters = []

    simulation.reporters.append(
        app.DCDReporter(f"traj_{i}.dcd", report_interval)
    )

    simulation.reporters.append(
        app.StateDataReporter(
            f"md_{i}.log",
            report_interval,
            step=True,
            potentialEnergy=True,
            temperature=True,
        )
    )

    simulation.step(steps_per_conformer)

print("All conformer MD runs finished.")