# Robot kinematics simulation

A MATLAB coursework project exploring task-space trajectories, inverse kinematics and forward kinematics for a three-degree-of-freedom robot model.

## What the script does

`three_dof_planar.m` defines link lengths and a desired end-effector trajectory, computes joint-space quantities, reconstructs Cartesian positions and visualizes the motion. The default script uses three 10 cm links, a 10-second motion and a 1 ms sampling interval.

```mermaid
flowchart LR
    A[Desired Cartesian trajectory] --> B[Inverse kinematics]
    B --> C[Joint positions and velocities]
    C --> D[Forward kinematics]
    D --> E[Plots and motion visualization]
```

## Run

Open MATLAB, set this repository as the current folder, and run:

```matlab
three_dof_planar
```

Inspect the parameters at the top of the script before changing the trajectory. Reachability and singular configurations matter: arbitrary endpoints can produce invalid inverse-kinematics expressions. MATLAB compatibility has not been revalidated across releases.

## Files

- `three_dof_planar.m`: main simulation.
- `backup.m`: alternate saved script.
- `θεωρητική_ανάλυση.pdf`: theoretical analysis in Greek.
- `εργαστηριακό_μέρος.pdf`: laboratory report in Greek.

This is a kinematic simulation, not a controller validated on physical hardware. Consult the accompanying reports for model assumptions and derivations.
