# coning_and_sculling

`coning_and_sculling` is a small Rust library for accumulating high-rate IMU
angular velocity and acceleration samples into lower-rate velocity and rotation
increments with coning and sculling corrections.

The crate exposes one public module, `coning_and_sculling`, containing the
`ConingAndSculling` accumulator. It uses `nalgebra::Vector3<f32>` for returned
vectors.

## Installation

Add this crate as a dependency from this repository:

```toml
[dependencies]
coning_and_sculling = { path = "../coning_sculling" }
```

The crate currently depends on:

```toml
nalgebra = "0.34"
```

## Usage

```rust
use coning_and_sculling::coning_and_sculling::ConingAndSculling;
use std::time::{Duration, Instant};

let start = Instant::now();
let mut corrector = ConingAndSculling::new(4, start);

for i in 1..=8 {
    let time = start + Duration::from_millis(i * 5);
    let angular_velocity = [0.0, 0.0, 0.1]; // rad/s, if using SI units
    let acceleration = [0.0, 0.0, 9.81];    // m/s^2, if using SI units

    if let Some((velocity_increment, rotation_vector)) =
        corrector.update(time, angular_velocity, acceleration)
    {
        println!("velocity increment: {velocity_increment:?}");
        println!("rotation vector: {rotation_vector:?}");
    }
}
```

`update` returns `None` while samples are still being accumulated. Once
`decimation_factor` samples have been processed, it returns:

- `velocity_increment`: integrated acceleration with rotational and sculling
  correction terms.
- `rotation_vector`: integrated angular velocity with coning correction terms.

Use consistent units for all samples. With angular velocity in radians per
second and acceleration in meters per second squared, the returned rotation
vector is in radians and the returned velocity increment is in meters per
second.

## API Overview

```rust
pub struct ConingAndSculling {
    pub decimation_factor: u32,
    pub sample: u32,
    pub time_prev: std::time::Instant,
    pub alpha: nalgebra::Vector3<f32>,
    pub delta_alpha: nalgebra::Vector3<f32>,
    pub nu: nalgebra::Vector3<f32>,
    pub delta_nu: nalgebra::Vector3<f32>,
    pub beta: nalgebra::Vector3<f32>,
    pub vel_scul: nalgebra::Vector3<f32>,
}
```

### `ConingAndSculling::new`

```rust
pub fn new(decimation_factor: u32, time: std::time::Instant) -> ConingAndSculling
```

Creates a new accumulator. Pass a positive `decimation_factor` to choose how
many IMU samples are combined into each output increment.

### `ConingAndSculling::update`

```rust
pub fn update(
    &mut self,
    time: std::time::Instant,
    angular_velocity: [f32; 3],
    acceleration: [f32; 3],
) -> Option<(nalgebra::Vector3<f32>, nalgebra::Vector3<f32>)>
```

Adds one IMU sample. The time step is computed from the difference between the
provided `time` and the previous sample time, so timestamps should be
monotonic.

### `ConingAndSculling::reset`

```rust
pub fn reset(&mut self, time: std::time::Instant)
```

Restarts accumulation from a new timestamp.

## Development

Run the test suite with:

```sh
cargo test
```

Run a quick compile check with:

```sh
cargo check
```

## Project Structure

```text
src/lib.rs
src/coning_and_sculling.rs
Cargo.toml
README.md
```

## License

No license file is currently included in this repository.
