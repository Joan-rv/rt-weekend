use crate::rt::{
    BvhNode, Camera, CameraDesc, CheckerTexture, ColorTexture, Dielectric, Hittable, Lambertian,
    Material, Metal, Sphere, Texture, render_pixel,
};
use glam::{U8Vec3, UVec2, Vec3};
use std::sync::Arc;
use wasm_bindgen::prelude::*;

#[wasm_bindgen]
#[derive(Debug, Clone, Copy)]
pub struct Rgb {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

impl From<Vec3> for Rgb {
    fn from(value: Vec3) -> Self {
        let U8Vec3 { x, y, z } = (value
            .map(|v| if v > 0.0 { v.sqrt() } else { v })
            .clamp(Vec3::splat(0.0), Vec3::splat(0.999))
            * 256.0)
            .as_u8vec3();
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

    pub fn add(&self, other: &WasmVec3) -> WasmVec3 {
        WasmVec3 {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }

    pub fn sub(&self, other: &WasmVec3) -> WasmVec3 {
        WasmVec3 {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }

    pub fn mul(&self, s: f32) -> WasmVec3 {
        WasmVec3 {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }

    #[wasm_bindgen(js_name = "mulVec")]
    pub fn mul_vec(&self, other: &WasmVec3) -> WasmVec3 {
        WasmVec3 {
            x: self.x * other.x,
            y: self.y * other.y,
            z: self.z * other.z,
        }
    }

    pub fn length(&self) -> f32 {
        (self.x * self.x + self.y * self.y + self.z * self.z).sqrt()
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
#[derive(Clone)]
pub struct TextureHandle(Arc<dyn Texture>);

#[wasm_bindgen]
#[derive(Clone)]
pub struct MaterialHandle(Arc<dyn Material>);

#[wasm_bindgen]
#[derive(Clone)]
pub struct HittableHandle(Arc<dyn Hittable>);

#[wasm_bindgen]
pub struct SceneBuilder;

#[wasm_bindgen]
impl SceneBuilder {
    #[wasm_bindgen(js_name = "createTextureSolidColor")]
    pub fn create_texture_solid_color(color: WasmVec3) -> TextureHandle {
        TextureHandle(Arc::new(ColorTexture {
            albedo: color.into(),
        }))
    }

    #[wasm_bindgen(js_name = "createTextureChecker")]
    pub fn create_texture_checker(
        scale: f32,
        even: &TextureHandle,
        odd: &TextureHandle,
    ) -> TextureHandle {
        TextureHandle(Arc::new(CheckerTexture {
            scale,
            even: even.0.clone(),
            odd: odd.0.clone(),
        }))
    }

    #[wasm_bindgen(js_name = "createMaterialLambertian")]
    pub fn create_material_lambertian(texture: &TextureHandle) -> MaterialHandle {
        MaterialHandle(Arc::new(Lambertian {
            texture: texture.0.clone(),
        }))
    }

    #[wasm_bindgen(js_name = "createMaterialMetal")]
    pub fn create_material_metal(albedo: WasmVec3, fuzz: f32) -> MaterialHandle {
        MaterialHandle(Arc::new(Metal {
            albedo: albedo.into(),
            fuzziness: fuzz,
        }))
    }

    #[wasm_bindgen(js_name = "createMaterialDielectric")]
    pub fn create_material_dielectric(refraction_index: f32) -> MaterialHandle {
        MaterialHandle(Arc::new(Dielectric { refraction_index }))
    }

    #[wasm_bindgen(js_name = "createSphere")]
    pub fn create_sphere(
        center: WasmVec3,
        radius: f32,
        material: &MaterialHandle,
    ) -> HittableHandle {
        HittableHandle(Arc::new(Sphere::stationary(
            center.into(),
            radius,
            material.0.clone(),
        )))
    }

    pub fn build(camera_desc: &WasmCameraDesc, hittables: Box<[HittableHandle]>) -> Scene {
        let mut world_list: Vec<Arc<dyn Hittable>> = Vec::new();
        for h in hittables.iter() {
            world_list.push(h.0.clone());
        }

        Scene {
            world: Box::new(BvhNode::create(&mut world_list)),
            camera: Camera::new(&camera_desc.clone().into()),
        }
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

