const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");
let image = ctx.createImageData(canvas.width, canvas.height);
let status = "stopped";

function draw() {
  if (status === "running") requestAnimationFrame(draw);
  ctx.putImageData(image, 0, 0);
}

function render() {
  image = ctx.createImageData(canvas.width, canvas.height);
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
      if (status === "running") {
        worker.postMessage([
          renderPos,
          Math.min(renderPos + batchSize, canvas.width * canvas.height),
        ]);
        renderPos += batchSize;
        if (renderPos >= canvas.width * canvas.height) {
          status = "stopped";
        }
      } else {
        worker.terminate();
      }
    };
  }
}

function buttonPress() {
  const button = document.getElementById("button");
  if (status === "running") {
    status = "stopped";
    button.innerHTML = "Render";
  } else {
    status = "running";
    button.innerHTML = "Stop";
    render();
  }
}
