import init, { SceneBuilder, CameraDesc, Vec3 } from "./pkg/rt_weekend.js";

async function run() {
  await init();

  const cameraDesc = new CameraDesc(
    100,
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
    const [x, y] = e.data;
    const pixel = scene.renderPixel(x, y);
    const { r, g, b } = pixel;
    pixel.free();
    postMessage({ r, g, b });
  };

  postMessage(null);
}

run();
