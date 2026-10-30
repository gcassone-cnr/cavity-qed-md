# Method and coordinate conventions

The recovered client evaluates a single-mode fixed-displacement QED-HF electronic surface using eT. i-PI propagates or samples the nuclei and one explicit photon coordinate. Electronic relaxation is carried out at each requested nuclear/photon configuration.

## Nuclear and photon coordinates

The i-PI system consists of all real nuclei followed by a pseudo-particle labelled `L`. Its x coordinate is the photon coordinate q; its mass is one in electron-mass atomic units. Its y/z coordinates and momenta stay zero. `L` is not sent to eT as a nucleus.

In standalone XYZ input to the workstation workflow, nuclear coordinates are in angstrom, but `L.x` is q in atomic units. Generated embedded i-PI bead arrays and the position trajectories marked `atomic_unit` use atomic units. Do not convert the photon coordinate using a nuclear angstrom-to-bohr factor. The supplied monomer PIMD initialization has q = 0, avoiding a nonzero mixed-unit initialization.

The client maps q to eT's fixed coherent displacement through

```math
z = q\sqrt{\omega/2}.
```

For the recovered sign convention the photon derivative is

```math
F_q = -\omega^2q - \omega\lambda\langle\boldsymbol\mu\cdot\mathbf e\rangle.
```

The electronic calculation supplies the complete energy for this implementation, including the photon harmonic and dipole-self-energy contributions. The driver must not add a second photon harmonic energy. The included fixed-coherent-state patch supplies the additional nuclear derivative terms on this path. These statements describe the supplied implementation; the force-check suite tests energy–derivative consistency on the target executable.

Use ordinary i-PI `ffsocket`. Adding a separate `ffcavphsocket` model would introduce additional cavity terms. This snapshot is designed for the supplied nonperiodic finite systems (`pbc=false`); it is not a validated periodic bulk-water implementation.

## Defaults differ between preserved workflows

| Setting | Workstation `workflow/` | `examples/monomer-pimd/` |
|---|---|---|
| Polarization | Laboratory z | Laboratory y |
| Frequency | 0.0187505 a.u. | 0.01666 a.u. |
| Coupling | 0.01 a.u. | 0.005 a.u. |
| Default basis | cc-pVDZ, tight SCF | aug-cc-pVDZ |
| Main purpose | Monomer/cluster NVT/NVE and one-bead Eckart NVE | Original monomer NVE and joint PIMD |

The two source trees are intentionally preserved separately. Their default outputs are not a controlled classical/quantum comparison. Match the physical parameters and electronic model before estimating a quantum correction.

## Statistical mechanics and interpretation

In the supplied PIMD input both nuclei and q are replicated over the ring polymer. A one-bead simulation and a many-bead simulation therefore change the treatment of both sectors. Separate nuclear, photonic and mixed quantum effects require additional matched protocols; they cannot be inferred from that comparison alone.

The eight-bead, 10 fs input is a software demonstration. Production equilibrium observables require bead, electronic-model and sampling convergence. A thermostatted PIMD trajectory is not an exact real-time quantum dynamics trajectory. The supplied IR and q analysis scripts require suitable real-time trajectories and correct sample spacing.

Global Eckart NVE constrains the nuclear center of mass and overall rotation relative to a reference, leaving internal molecular/cluster motion active and excluding the photon from the nuclear group. With fixed cavity polarization, this changes the allowed rotational motion. The recovered Eckart path is restricted to one bead and NVE; Eckart PIMD/NVT are not advertised here.

Charged clusters, multiple modes, loss models, driven cavities, PLUMED integration and quantitative reaction-rate calculations need their own explicit source, conventions and validation before being added to a release.
