extern crate nalgebra as na;
use na::Vector3;
use std::time;

pub struct ConingAndSculling {
    pub decimation_factor: u32,
    pub sample: u32,
    pub time_prev: time::Instant,
    pub alpha: Vector3<f32>,
    pub delta_alpha: Vector3<f32>,
    pub nu: Vector3<f32>,
    pub delta_nu: Vector3<f32>,
    pub beta: Vector3<f32>,
    pub vel_scul: Vector3<f32>,
}

impl ConingAndSculling {
    pub fn new(decimation_factor: u32, time: time::Instant) -> ConingAndSculling {
        assert!(
            decimation_factor > 0,
            "decimation_factor must be greater than zero"
        );

        Self {
            decimation_factor,
            sample: 1,
            time_prev: time,
            alpha: Vector3::zeros(),
            delta_alpha: Vector3::zeros(),
            nu: Vector3::zeros(),
            delta_nu: Vector3::zeros(),
            beta: Vector3::zeros(),
            vel_scul: Vector3::zeros(),
        }
    }

    pub fn reset(&mut self, time: time::Instant) {
        self.time_prev = time;
        self.sample = 1;
        self.alpha = Vector3::zeros();
        self.delta_alpha = Vector3::zeros();
        self.nu = Vector3::zeros();
        self.delta_nu = Vector3::zeros();
        self.beta = Vector3::zeros();
        self.vel_scul = Vector3::zeros();
    }

    pub fn update(
        &mut self,
        time: time::Instant,
        angular_velocity: [f32; 3],
        acceleration: [f32; 3],
    ) -> Option<(Vector3<f32>, Vector3<f32>)> {
        assert!(
            self.decimation_factor > 0,
            "decimation_factor must be greater than zero"
        );

        let delta_time: f32 = (time - self.time_prev).as_secs_f32();
        // println!("delta_time: {} | sample: {} | ", delta_time, self.sample);

        if self.sample <= self.decimation_factor {
            self.time_prev = time;
            let alpha_prev: Vector3<f32> = self.alpha;
            let nu_prev: Vector3<f32> = self.nu;
            let delta_alpha_prev: Vector3<f32> = self.delta_alpha;
            let delta_nu_prev: Vector3<f32> = self.delta_nu;
            let beta_prev: Vector3<f32> = self.beta;

            self.delta_alpha = Vector3::from(angular_velocity) * delta_time;
            self.delta_nu = Vector3::from(acceleration) * delta_time;

            self.alpha = alpha_prev + self.delta_alpha;
            self.nu = nu_prev + self.delta_nu;

            let delta_beta: Vector3<f32> =
                0.5 * (alpha_prev + delta_alpha_prev / 6.0).cross(&self.delta_alpha);
            self.beta = beta_prev + delta_beta;

            let vel_scul_prev: Vector3<f32> = self.vel_scul;
            let mut first_factor: Vector3<f32> = alpha_prev.cross(&(delta_alpha_prev / 6.0));
            first_factor = first_factor.cross(&self.delta_nu);

            let second_factor: Vector3<f32> =
                (nu_prev + delta_nu_prev / 6.0).cross(&self.delta_alpha);
            let delta_vel_scul: Vector3<f32> = 0.5 * (first_factor + second_factor);
            self.vel_scul = vel_scul_prev + delta_vel_scul;

            self.sample += 1;
        }

        if self.sample > self.decimation_factor {
            self.sample = 1;
            let vel_rot: Vector3<f32> = 0.5 * self.alpha.cross(&self.nu);
            let vel_imu: Vector3<f32> = self.nu + vel_rot + self.vel_scul;
            let rot_vec_imu: Vector3<f32> = self.alpha + self.beta;
            self.reset(time);
            Some((vel_imu, rot_vec_imu))
        } else {
            None
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;

    fn assert_vec_near(actual: Vector3<f32>, expected: Vector3<f32>) {
        let error = (actual - expected).norm();
        assert!(
            error < 1e-6,
            "expected {expected:?}, got {actual:?}, error {error}"
        );
    }

    #[test]
    fn includes_coning_correction_from_current_sample_increment() {
        let t_0 = time::Instant::now();
        let mut coning_and_sculling = ConingAndSculling::new(2, t_0);

        assert_eq!(
            coning_and_sculling.update(
                t_0 + time::Duration::from_secs(1),
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
            ),
            None
        );

        let (_, rot_vec_imu) = coning_and_sculling
            .update(
                t_0 + time::Duration::from_secs(2),
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
            )
            .expect("second sample should complete the decimation window");

        assert_vec_near(rot_vec_imu, Vector3::new(1.0, 1.0, 7.0 / 12.0));
    }

    #[test]
    fn reset_clears_sculling_correction_between_windows() {
        let t_0 = time::Instant::now();
        let mut coning_and_sculling = ConingAndSculling::new(2, t_0);

        assert_eq!(
            coning_and_sculling.update(
                t_0 + time::Duration::from_secs(1),
                [1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
            ),
            None
        );

        let (vel_imu, _) = coning_and_sculling
            .update(
                t_0 + time::Duration::from_secs(2),
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
            )
            .expect("second sample should complete the first decimation window");
        assert!(vel_imu.z != 0.0, "test setup should generate sculling");

        assert_eq!(
            coning_and_sculling.update(
                t_0 + time::Duration::from_secs(3),
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
            ),
            None
        );

        let (vel_imu, rot_vec_imu) = coning_and_sculling
            .update(
                t_0 + time::Duration::from_secs(4),
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
            )
            .expect("second sample should complete the second decimation window");

        assert_vec_near(vel_imu, Vector3::zeros());
        assert_vec_near(rot_vec_imu, Vector3::zeros());
    }

    #[test]
    #[should_panic(expected = "decimation_factor must be greater than zero")]
    fn rejects_zero_decimation_factor() {
        let _ = ConingAndSculling::new(0, time::Instant::now());
    }

    #[test]
    #[should_panic(expected = "decimation_factor must be greater than zero")]
    fn rejects_externally_mutated_zero_decimation_factor() {
        let t_0 = time::Instant::now();
        let mut coning_and_sculling = ConingAndSculling::new(2, t_0);

        coning_and_sculling.decimation_factor = 0;

        let _ = coning_and_sculling.update(
            t_0 + time::Duration::from_secs(1),
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        );
    }

    #[test]
    fn public_fields_remain_accessible() {
        let t_0 = time::Instant::now();
        let mut coning_and_sculling = ConingAndSculling::new(2, t_0);

        coning_and_sculling.sample = 1;
        coning_and_sculling.time_prev = t_0;
        coning_and_sculling.alpha = Vector3::zeros();
        coning_and_sculling.delta_alpha = Vector3::zeros();
        coning_and_sculling.nu = Vector3::zeros();
        coning_and_sculling.delta_nu = Vector3::zeros();
        coning_and_sculling.beta = Vector3::zeros();
        coning_and_sculling.vel_scul = Vector3::zeros();

        assert_eq!(coning_and_sculling.decimation_factor, 2);
    }
}
