import init, { SceneBuilder, CameraDesc, Vec3 } from "./pkg/rt_weekend.js";

async function run() {
  await init();

  const cameraDesc = new CameraDesc(
    800,
    1.0,
    new Vec3(13.0, 2.0, 3.0),
    new Vec3(0.0, 0.0, 0.0),
    (0.6 * Math.PI) / 180.0,
    10.0,
    (20.0 * Math.PI) / 180.0,
    100,
    10,
  );
  const sceneBuilder = new SceneBuilder(cameraDesc);
  sceneBuilder.addSphereStationary(new Vec3(0.0, 0.0, 0.0), 1.0);
  const scene = sceneBuilder.build();

  onmessage = (e) => {
    let pixels = [];
    const [start, end] = e.data;
    for (let i = start; i < end; ++i) {
      const x = i % 800;
      const y = Math.floor(i / 800);
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
