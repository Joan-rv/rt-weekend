const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");
let image = ctx.createImageData(canvas.width, canvas.height);
let done;

function draw() {
  ctx.putImageData(image, 0, 0);
  if (!done) requestAnimationFrame(draw);
}

async function run() {
  done = false;
  requestAnimationFrame(draw);

  let renderPos = 0;
  const batchSize = 10;
  for (let i = 0; i < navigator.hardwareConcurrency; ++i) {
    const worker = new Worker("worker.js", { type: "module" });
    worker.onmessage = (e) => {
      if (e.data) {
        for (let {
          coords: { x, y },
          rgb: { r, g, b },
        } of e.data) {
          image.data[4 * (y * canvas.width + x)] = r;
          image.data[4 * (y * canvas.width + x) + 1] = g;
          image.data[4 * (y * canvas.width + x) + 2] = b;
          image.data[4 * (y * canvas.width + x) + 3] = 0xff;
        }
      }
      if (!done) {
        worker.postMessage([
          renderPos,
          Math.min(renderPos + batchSize, canvas.width * canvas.height),
        ]);
        renderPos += batchSize;
        if (renderPos >= canvas.width * canvas.height) {
          done = true;
        }
      } else {
        worker.terminate();
      }
    };
  }
}

run();
