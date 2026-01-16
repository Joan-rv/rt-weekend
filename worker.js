import init, { SceneBuilder, CameraDesc, Vec3 } from "./pkg/rt_weekend.js";

let seed = 123456789;
function random() {
  seed = (Math.imul(1664525, seed) + 1013904223) | 0;
  return (seed >>> 0) / 4294967296;
}

function randomRange(min, max) {
  return min + (max - min) * random();
}

function randomVec3() {
  return new Vec3(random(), random(), random());
}

function randomVec3Range(min, max) {
  return new Vec3(
    randomRange(min, max),
    randomRange(min, max),
    randomRange(min, max),
  );
}

async function run() {
  await init();

  seed = 123456789;

  const imageWidth = 800;
  const aspectRatio = 1.0;

  const hittables = [];

  const checkerTexture = SceneBuilder.createTextureChecker(
    0.32,
    SceneBuilder.createTextureSolidColor(new Vec3(0.2, 0.3, 0.1)),
    SceneBuilder.createTextureSolidColor(new Vec3(0.9, 0.9, 0.9)),
  );
  const groundMat = SceneBuilder.createMaterialLambertian(checkerTexture);
  hittables.push(
    SceneBuilder.createSphere(new Vec3(0.0, -1000.0, 0.0), 1000.0, groundMat),
  );

  for (let a = -11; a < 11; a++) {
    for (let b = -11; b < 11; b++) {
      const choose_mat = random();
      const center = new Vec3(a + 0.9 * random(), 0.2, b + 0.9 * random());

      if (center.sub(new Vec3(4.0, 0.2, 0.0)).length() > 0.9) {
        if (choose_mat < 0.8) {
          const albedo = randomVec3().mulVec(randomVec3());
          const sphereMat = SceneBuilder.createMaterialLambertian(
            SceneBuilder.createTextureSolidColor(albedo),
          );
          hittables.push(SceneBuilder.createSphere(center, 0.2, sphereMat));
        } else if (choose_mat < 0.95) {
          const albedo = randomVec3Range(0.5, 1.0);
          const fuzz = randomRange(0.0, 0.5);
          const sphereMat = SceneBuilder.createMaterialMetal(albedo, fuzz);
          hittables.push(SceneBuilder.createSphere(center, 0.2, sphereMat));
        } else {
          const sphereMat = SceneBuilder.createMaterialDielectric(1.5);
          hittables.push(SceneBuilder.createSphere(center, 0.2, sphereMat));
        }
      }
    }
  }

  const mat1 = SceneBuilder.createMaterialDielectric(1.5);
  hittables.push(SceneBuilder.createSphere(new Vec3(0.0, 1.0, 0.0), 1.0, mat1));

  const mat2 = SceneBuilder.createMaterialLambertian(
    SceneBuilder.createTextureSolidColor(new Vec3(0.4, 0.2, 0.1)),
  );
  hittables.push(
    SceneBuilder.createSphere(new Vec3(-4.0, 1.0, 0.0), 1.0, mat2),
  );

  const mat3 = SceneBuilder.createMaterialMetal(new Vec3(0.7, 0.6, 0.5), 0.0);
  hittables.push(SceneBuilder.createSphere(new Vec3(4.0, 1.0, 0.0), 1.0, mat3));

  const cameraDesc = {
    imageWidth,
    aspectRatio,
    lookFrom: { x: 13.0, y: 2.0, z: 3.0 },
    lookAt: { x: 0.0, y: 0.0, z: 0.0 },
    defocusAngle: (0.6 * Math.PI) / 180.0,
    focusDist: 10.0,
    fovy: (20.0 * Math.PI) / 180.0,
    samplesPerPx: 500,
    maxDepth: 50,
  };
  const scene = SceneBuilder.build(cameraDesc, hittables);

  onmessage = (e) => {
    let pixels = [];
    const [start, end] = e.data;
    for (let i = start; i < end; ++i) {
      const x = i % imageWidth;
      const y = Math.floor(i / imageWidth);
      const pixel = scene.renderPixel(x, y);
      const { r, g, b } = pixel;
      pixel.free();
      pixels.push({ coords: { x, y }, rgb: { r, g, b } });
    }
    postMessage(pixels);
  };

  postMessage(null);
}

run();
