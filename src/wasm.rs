use crate::rt::{BvhNode, Camera, CameraDesc, Hittable, Lambertian, Sphere, render_pixel};
use glam::{UVec2, Vec3};
use std::sync::Arc;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
#[derive(Debug, Clone, Copy)]
pub struct Rgb {
    pub r: f32,
    pub g: f32,
    pub b: f32,
}

impl From<Vec3> for Rgb {
    fn from(value: Vec3) -> Self {
        let Vec3 { x, y, z } = value;
        Self { r: x, g: y, b: z }
    }
}

#[wasm_bindgen(js_name = "Vec3")]
#[derive(Debug, Clone, Copy)]
pub struct WasmVec3 {
    pub x: f32,
    pub y: f32,
    pub z: f32,
}

#[wasm_bindgen(js_class = "Vec3")]
impl WasmVec3 {
    #[wasm_bindgen(constructor)]
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
    }
}

impl From<Vec3> for WasmVec3 {
    fn from(value: Vec3) -> Self {
        let Vec3 { x, y, z } = value;
        Self { x, y, z }
    }
}

impl From<WasmVec3> for Vec3 {
    fn from(value: WasmVec3) -> Self {
        let WasmVec3 { x, y, z } = value;
        Self { x, y, z }
    }
}

#[wasm_bindgen]
pub struct Scene {
    world: Box<dyn Hittable>,
    camera: Camera,
}

#[wasm_bindgen]
impl Scene {
    #[wasm_bindgen(js_name = "renderPixel")]
    pub fn render_pixel(&self, x: u32, y: u32) -> Rgb {
        render_pixel(&self.camera, &*self.world, UVec2::new(x, y)).into()
    }
}

#[wasm_bindgen]
pub struct SceneBuilder {
    hittables: Vec<Arc<dyn Hittable>>,
    camera: Camera,
}

#[wasm_bindgen]
impl SceneBuilder {
    #[wasm_bindgen(constructor)]
    pub fn new(desc: &WasmCameraDesc) -> Self {
        console_error_panic_hook::set_once();
        Self {
            hittables: Vec::new(),
            camera: Camera::new(&desc.clone().into()),
        }
    }

    #[wasm_bindgen]
    pub fn build(mut self) -> Scene {
        Scene {
            world: Box::new(BvhNode::create(&mut self.hittables)),
            camera: self.camera,
        }
    }

    #[wasm_bindgen(js_name = "addSphereStationary")]
    pub fn add_sphere_stationary(&mut self, center: WasmVec3, radius: f32) {
        self.hittables.push(Arc::new(Sphere::stationary(
            center.into(),
            radius,
            Arc::new(Lambertian::from_color(Vec3::new(1.0, 0.5, 0.0))),
        )));
    }
}

#[wasm_bindgen(js_name = "CameraDesc")]
#[derive(Debug, Clone)]
pub struct WasmCameraDesc {
    pub image_width: u32,
    pub aspect_ratio: f32,

    pub look_from: WasmVec3,
    pub look_at: WasmVec3,

    pub defocus_angle: f32,
    pub focus_dist: f32,
    pub fovy: f32,

    pub samples_per_px: u32,
    pub max_depth: u32,
}

#[wasm_bindgen(js_class = "CameraDesc")]
impl WasmCameraDesc {
    #[wasm_bindgen(constructor)]
    pub fn new(
        image_width: u32,
        aspect_ratio: f32,
        look_from: WasmVec3,
        look_at: WasmVec3,
        defocus_angle: f32,
        focus_dist: f32,
        fovy: f32,
        samples_per_px: u32,
        max_depth: u32,
    ) -> Self {
        Self {
            image_width,
            aspect_ratio,
            look_from,
            look_at,
            defocus_angle,
            focus_dist,
            fovy,
            samples_per_px,
            max_depth,
        }
    }
}

impl From<WasmCameraDesc> for CameraDesc {
    fn from(w: WasmCameraDesc) -> Self {
        Self {
            image_width: w.image_width,
            aspect_ratio: w.aspect_ratio,
            look_from: w.look_from.into(),
            look_at: w.look_at.into(),
            defocus_angle: w.defocus_angle,
            focus_dist: w.focus_dist,
            fovy: w.fovy,
            samples_per_px: w.samples_per_px,
            max_depth: w.max_depth,
        }
    }
}
