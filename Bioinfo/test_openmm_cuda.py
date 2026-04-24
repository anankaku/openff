import openmm

print("OpenMM version:", openmm.version.version)

print("Available OpenMM platforms:")
for i in range(openmm.Platform.getNumPlatforms()):
    p = openmm.Platform.getPlatform(i)
    print("-", p.getName())

try:
    cuda = openmm.Platform.getPlatformByName("CUDA")
    print("CUDA platform available")
except Exception as e:
    print("CUDA platform NOT available")
    print(e)