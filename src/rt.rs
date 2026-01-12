use crate::rng::fast_random;
use glam::{UVec2, Vec2, Vec3};
use image::Rgb32FImage;
use std::f32::{self, consts::PI};
use std::{ops::Index, sync::Arc};

pub const INFINITY: f32 = std::f32::MAX;

#[derive(Debug, Clone, Copy)]
pub struct Interval {
    pub min: f32,
    pub max: f32,
}

impl Interval {
    pub const EMPTY: Self = Self {
        min: INFINITY,
        max: -INFINITY,
    };
    pub const UNIVERSE: Self = Self {
        min: -INFINITY,
        max: INFINITY,
    };

    pub fn new(min: f32, max: f32) -> Self {
        Self { min, max }
    }

    pub fn from_intervals(a: Self, b: Self) -> Self {
        let min = a.min.min(b.min);
        let max = a.max.max(b.max);
        Self { min, max }
    }

    pub fn size(&self) -> f32 {
        self.max - self.min
    }

    pub fn contains(&self, x: f32) -> bool {
        self.min <= x && x <= self.max
    }
    pub fn surrounds(&self, x: f32) -> bool {
        self.min < x && x < self.max
    }

    pub fn expand(&self, delta: f32) -> Self {
        let padding = delta / 2.0;
        Interval {
            min: self.min - padding,
            max: self.max + padding,
        }
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Aabb {
    pub x: Interval,
    pub y: Interval,
    pub z: Interval,
}

impl Aabb {
    pub const EMPTY: Self = Self::new(Interval::EMPTY, Interval::EMPTY, Interval::EMPTY);

    pub const fn new(x: Interval, y: Interval, z: Interval) -> Self {
        Self { x, y, z }
    }

    pub fn from_points(a: Vec3, b: Vec3) -> Self {
        let x = if a[0] < b[0] {
            Interval::new(a[0], b[0])
        } else {
            Interval::new(b[0], a[0])
        };
        let y = if a[1] < b[1] {
            Interval::new(a[1], b[1])
        } else {
            Interval::new(b[1], a[1])
        };
        let z = if a[2] < b[2] {
            Interval::new(a[2], b[2])
        } else {
            Interval::new(b[2], a[2])
        };
        Self { x, y, z }
    }

    pub fn from_aabbs(a: Self, b: Self) -> Self {
        Self {
            x: Interval::from_intervals(a.x, b.x),
            y: Interval::from_intervals(a.y, b.y),
            z: Interval::from_intervals(a.z, b.z),
        }
    }

    pub fn longest_axis(&self) -> usize {
        if self.x.size() > self.y.size() {
            if self.x.size() > self.z.size() { 0 } else { 2 }
        } else if self.y.size() > self.z.size() {
            1
        } else {
            0
        }
    }

    pub fn hit(&self, ray: &Ray, mut valid_t: Interval) -> bool {
        for axis in 0..3 {
            let interval = [&self.x, &self.y, &self.z][axis];
            let dir_inv = 1.0 / ray.direction[axis];
            let t0 = (interval.min - ray.origin[axis]) * dir_inv;
            let t1 = (interval.max - ray.origin[axis]) * dir_inv;
            if t0 < t1 {
                valid_t.min = valid_t.min.max(t0);
                valid_t.max = valid_t.max.min(t1);
            } else {
                valid_t.min = valid_t.min.max(t1);
                valid_t.max = valid_t.max.min(t0);
            }

            if valid_t.max <= valid_t.min {
                return false;
            }
        }
        true
    }
}

impl Index<usize> for Aabb {
    type Output = Interval;

    fn index(&self, index: usize) -> &Self::Output {
        match index {
            0 => &self.x,
            1 => &self.y,
            2 => &self.z,
            _ => panic!("{index} is out of range for Aabb"),
        }
    }
}

pub fn random_unit_vec3() -> Vec3 {
    loop {
        let v = 2.0 * fast_random::<Vec3>() - 1.0;
        let len_sq = v.length_squared();
        if len_sq < 1.0
            && let Some(v) = v.try_normalize()
        {
            return v;
        }
    }
}

pub fn random_unit_disk_vec3() -> Vec3 {
    loop {
        let v = 2.0 * fast_random::<Vec3>() - 1.0;
        if v.length_squared() < 1.0 {
            return v;
        }
    }
}

#[derive(Debug, Clone)]
pub struct Ray {
    pub origin: Vec3,
    pub direction: Vec3,
    pub time: f32,
}

impl Ray {
    pub fn at(&self, t: f32) -> Vec3 {
        self.origin + self.direction * t
    }
}

#[derive(Clone)]
pub struct HitRecord<'material> {
    pub point: Vec3,
    pub normal: Vec3,
    pub tex_coords: Vec2,
    pub material: &'material dyn Material,
    pub t: f32,
    pub front_face: bool,
}

pub trait Hittable: Sync + Send {
    fn hit(&self, ray: &Ray, valid_t: Interval) -> Option<HitRecord<'_>>;
    fn bounding_box(&self) -> Aabb;
}

pub trait Material: Sync + Send {
    fn scatter(&self, ray: &Ray, record: &HitRecord) -> Option<(Vec3, Ray)>;
}

#[derive(Debug, Clone)]
pub struct Camera {
    image_width: u32,
    image_height: u32,

    center: Vec3,
    pixel00: Vec3,
    pixel_du: Vec3,
    pixel_dv: Vec3,

    defocus_angle: f32,
    defocus_disk_u: Vec3,
    defocus_disk_v: Vec3,

    samples_per_px: u32,
    max_depth: u32,
}

#[derive(Debug, Clone)]
pub struct CameraDesc {
    pub image_width: u32,
    pub aspect_ratio: f32,

    pub look_from: Vec3,
    pub look_at: Vec3,

    pub defocus_angle: f32,
    pub focus_dist: f32,

    pub fovy: f32,

    pub samples_per_px: u32,
    pub max_depth: u32,
}

impl Camera {
    pub fn new(desc: &CameraDesc) -> Self {
        let image_height = (desc.image_width as f32 / desc.aspect_ratio) as u32;
        let h = (desc.fovy / 2.0).tan();

        let viewport_height = 2.0 * h * desc.focus_dist;
        let viewport_width = viewport_height * desc.image_width as f32 / image_height as f32;

        let w = (desc.look_from - desc.look_at).normalize();
        let u = Vec3::Y.cross(w).normalize();
        let v = w.cross(u);

        let viewport_u = viewport_width * u;
        let viewport_v = viewport_height * -v;
        let pixel_du = viewport_u / (desc.image_width as f32);
        let pixel_dv = viewport_v / (image_height as f32);

        let viewport_top_left =
            desc.look_from - (desc.focus_dist * w) - viewport_u / 2.0 - viewport_v / 2.0;
        let pixel00 = viewport_top_left + 0.5 * (pixel_du + pixel_dv);

        let defocus_radius = desc.focus_dist * (desc.defocus_angle / 2.0).tan();
        let defocus_disk_u = u * defocus_radius;
        let defocus_disk_v = v * defocus_radius;
        Self {
            image_width: desc.image_width,
            image_height,
            center: desc.look_from,
            pixel00,
            pixel_du,
            pixel_dv,

            defocus_angle: desc.defocus_angle,
            defocus_disk_u,
            defocus_disk_v,

            samples_per_px: desc.samples_per_px,
            max_depth: desc.max_depth,
        }
    }

    pub fn width(&self) -> u32 {
        self.image_width
    }
    pub fn height(&self) -> u32 {
        self.image_height
    }

    fn get_ray(&self, coords: UVec2) -> Ray {
        let x = coords.x as f32 + fast_random::<f32>() - 0.5;
        let y = coords.y as f32 + fast_random::<f32>() - 0.5;
        let px_loc = self.pixel00 + x * self.pixel_du + y * self.pixel_dv;
        let origin = if self.defocus_angle > 0.0 {
            let offset = random_unit_disk_vec3();
            self.center + offset.x * self.defocus_disk_u + offset.y * self.defocus_disk_v
        } else {
            self.center
        };
        let time = fast_random();
        Ray {
            origin,
            direction: px_loc - origin,
            time,
        }
    }
}

fn ray_color(ray: &Ray, depth: u32, world: &dyn Hittable) -> Vec3 {
    if depth == 0 {
        return Vec3::ZERO;
    }

    let pos_interval = Interval {
        min: 0.001,
        max: INFINITY,
    };
    if let Some(record) = world.hit(ray, pos_interval) {
        if let Some((attenuation, scattered)) = record.material.scatter(ray, &record) {
            return attenuation * ray_color(&scattered, depth - 1, world);
        }
        return Vec3::X;
    }

    let dir = ray.direction.normalize();
    let a = 0.5 * (dir.y + 1.0);
    (1.0 - a) * Vec3::new(1.0, 1.0, 1.0) + a * Vec3::new(0.5, 0.7, 1.0)
}

pub fn render_pixel(camera: &Camera, world: &dyn Hittable, coords: UVec2) -> Vec3 {
    let mut color = Vec3::ZERO;
    for _ in 0..camera.samples_per_px {
        color += ray_color(&camera.get_ray(coords), camera.max_depth, world);
    }
    color / (camera.samples_per_px as f32)
}

#[derive(Clone)]
pub struct Sphere {
    pub center: Ray,
    pub radius: f32,
    pub material: Arc<dyn Material>,
}

impl Sphere {
    pub fn stationary(center: Vec3, radius: f32, material: Arc<dyn Material>) -> Self {
        Self {
            center: Ray {
                origin: center,
                direction: Vec3::ZERO,
                time: 0.0,
            },
            radius,
            material,
        }
    }

    pub fn moving(center1: Vec3, center2: Vec3, radius: f32, material: Arc<dyn Material>) -> Self {
        Self {
            center: Ray {
                origin: center1,
                direction: center2 - center1,
                time: 0.0,
            },
            radius,
            material,
        }
    }
}

impl Hittable for Sphere {
    fn hit(&self, ray: &Ray, valid_t: Interval) -> Option<HitRecord<'_>> {
        let center = self.center.at(ray.time);
        let oc = center - ray.origin;
        let a = ray.direction.length_squared();
        let h = ray.direction.dot(oc);
        let c = oc.length_squared() - self.radius * self.radius;
        let discriminant = h * h - a * c;
        if discriminant < 0.0 {
            return None;
        }

        let discriminant_sqrt = discriminant.sqrt();
        let mut root = (h - discriminant_sqrt) / a;
        if !valid_t.surrounds(root) {
            root = (h + discriminant_sqrt) / a;
            if !valid_t.surrounds(root) {
                return None;
            }
        }

        let point = ray.at(root);
        let outward_normal = (point - center) / self.radius;
        let front_face = ray.direction.dot(outward_normal) < 0.0;
        let normal = if front_face {
            outward_normal
        } else {
            -outward_normal
        };
        let theta = f32::acos(-outward_normal.y);
        let phi = f32::atan2(-outward_normal.z, outward_normal.x) + PI;
        let u = phi / (2.0 * PI);
        let v = theta / PI;
        let tex_coords = Vec2::new(u, v);

        Some(HitRecord {
            point,
            normal,
            tex_coords,
            t: root,
            material: &*self.material,
            front_face,
        })
    }

    fn bounding_box(&self) -> Aabb {
        let rvec = Vec3::splat(self.radius);
        let box0 = Aabb::from_points(self.center.at(0.0) - rvec, self.center.at(0.0) + rvec);
        let box1 = Aabb::from_points(self.center.at(1.0) - rvec, self.center.at(1.0) + rvec);
        Aabb::from_aabbs(box0, box1)
    }
}

impl Hittable for Vec<Box<dyn Hittable>> {
    fn hit(&self, ray: &Ray, valid_t: Interval) -> Option<HitRecord<'_>> {
        let mut closest_hit = None;
        let mut closest_interval = valid_t;
        for hittable in self {
            if let Some(record) = hittable.hit(ray, closest_interval) {
                closest_interval.max = record.t;
                closest_hit = Some(record);
            }
        }
        closest_hit
    }

    fn bounding_box(&self) -> Aabb {
        let mut bbox = Aabb::EMPTY;
        for hittable in self {
            bbox = Aabb::from_aabbs(bbox, hittable.bounding_box());
        }
        bbox
    }
}

type BvhLink = Arc<dyn Hittable>;
#[derive(Clone)]
pub struct BvhNode {
    left: BvhLink,
    right: BvhLink,
    bbox: Aabb,
}

impl BvhNode {
    pub fn create(hittables: &mut [Arc<dyn Hittable>]) -> Self {
        if hittables.len() == 0 {
            panic!("must create BvhNode with at least one hittable");
        }

        let mut bbox = Aabb::EMPTY;
        for hittable in &mut *hittables {
            bbox = Aabb::from_aabbs(bbox, hittable.bounding_box());
        }
        if hittables.len() == 1 {
            return BvhNode {
                left: hittables[0].clone(),
                right: hittables[0].clone(),
                bbox,
            };
        }
        if hittables.len() == 2 {
            return BvhNode {
                left: hittables[0].clone(),
                right: hittables[1].clone(),
                bbox,
            };
        }

        let axis = bbox.longest_axis();
        hittables.sort_unstable_by(|a, b| {
            a.bounding_box()[axis]
                .min
                .total_cmp(&b.bounding_box()[axis].min)
        });
        let mid = hittables.len() / 2;
        let (left, right) = hittables.split_at_mut(mid);
        let left = Arc::new(BvhNode::create(left));
        let right = Arc::new(BvhNode::create(right));
        BvhNode { left, right, bbox }
    }
}

impl Hittable for BvhNode {
    fn hit(&self, ray: &Ray, valid_t: Interval) -> Option<HitRecord<'_>> {
        if !self.bbox.hit(ray, valid_t) {
            return None;
        }

        let left = self.left.hit(ray, valid_t);
        let max_t = if let Some(ref record) = left {
            record.t
        } else {
            valid_t.max
        };
        let right = self.right.hit(ray, Interval::new(valid_t.min, max_t));

        right.or(left)
    }

    fn bounding_box(&self) -> Aabb {
        self.bbox
    }
}

#[derive(Clone)]
pub struct Lambertian {
    pub texture: Arc<dyn Texture>,
}

impl Lambertian {
    pub fn from_color(albedo: Vec3) -> Self {
        Self {
            texture: Arc::new(ColorTexture { albedo }),
        }
    }
}

impl Material for Lambertian {
    fn scatter(&self, ray: &Ray, record: &HitRecord) -> Option<(Vec3, Ray)> {
        let mut out_dir = random_unit_vec3() + record.normal;
        if out_dir.length_squared() < 1e-4 {
            out_dir = record.normal;
        }
        Some((
            self.texture.value(record.tex_coords, record.point),
            Ray {
                direction: out_dir,
                origin: record.point,
                time: ray.time,
            },
        ))
    }
}

#[derive(Debug, Clone)]
pub struct Metal {
    pub albedo: Vec3,
    pub fuzziness: f32,
}

impl Material for Metal {
    fn scatter(&self, ray: &Ray, record: &HitRecord) -> Option<(Vec3, Ray)> {
        Some((
            self.albedo,
            Ray {
                origin: record.point,
                direction: ray.direction.reflect(record.normal)
                    + self.fuzziness * random_unit_vec3(),
                time: ray.time,
            },
        ))
    }
}

#[derive(Debug, Clone)]
pub struct Dielectric {
    pub refraction_index: f32,
}

impl Dielectric {
    fn reflectance(&self, cos_theta: f32) -> f32 {
        // Schlick's approximation for reflectance
        let r0 = (1.0 - self.refraction_index) / (1.0 + self.refraction_index);
        let r0 = r0 * r0;
        r0 + (1.0 - r0) * (1.0 - cos_theta).powi(5)
    }
}

impl Material for Dielectric {
    fn scatter(&self, ray: &Ray, record: &HitRecord) -> Option<(Vec3, Ray)> {
        let refraction_index = if record.front_face {
            1.0 / self.refraction_index
        } else {
            self.refraction_index
        };
        let unit_direction = ray.direction.normalize();
        let cos_theta = (-unit_direction).dot(record.normal).min(1.0);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        let cannot_refract = refraction_index * sin_theta > 1.0;
        let direction = if cannot_refract || self.reflectance(cos_theta) > fast_random() {
            unit_direction.reflect(record.normal)
        } else {
            unit_direction.refract(record.normal, refraction_index)
        };

        Some((
            Vec3::ONE,
            Ray {
                origin: record.point,
                direction,
                time: ray.time,
            },
        ))
    }
}

pub trait Texture: Send + Sync {
    fn value(&self, tex_coords: Vec2, point: Vec3) -> Vec3;
}

#[derive(Debug, Clone)]
pub struct ColorTexture {
    pub albedo: Vec3,
}

impl Texture for ColorTexture {
    fn value(&self, _tex_coords: Vec2, _point: Vec3) -> Vec3 {
        self.albedo
    }
}

#[derive(Clone)]
pub struct CheckerTexture {
    pub scale: f32,
    pub even: Arc<dyn Texture>,
    pub odd: Arc<dyn Texture>,
}

impl Texture for CheckerTexture {
    fn value(&self, tex_coords: Vec2, point: Vec3) -> Vec3 {
        let point_int = (point / self.scale).floor().as_ivec3();
        let is_even = point_int.element_sum() % 2 == 0;
        if is_even {
            self.even.value(tex_coords, point)
        } else {
            self.odd.value(tex_coords, point)
        }
    }
}

impl Texture for Rgb32FImage {
    fn value(&self, tex_coords: Vec2, _point: Vec3) -> Vec3 {
        let u = tex_coords.x.clamp(0.0, 1.0);
        let v = 1.0 - tex_coords.y.clamp(0.0, 1.0);
        let x = (u * self.width() as f32) as u32;
        let y = (v * self.height() as f32) as u32;
        let px: Vec3 = self.get_pixel(x, y).0.into();
        px.powf(2.0)
    }
}
