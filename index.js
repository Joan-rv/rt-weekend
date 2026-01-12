const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");
let image = ctx.createImageData(100, 100);
let done = false;

function draw() {
  ctx.putImageData(image, 0, 0);
  if (!done) requestAnimationFrame(draw);
}

async function run() {
  requestAnimationFrame(draw);

  let x = 0;
  let y = 0;
  let workers = [];
  for (let i = 0; i < navigator.hardwareConcurrency; ++i) {
    const worker = new Worker("worker.js", { type: "module" });
    worker.onmessage = (e) => {
      if (e.data) {
        const {
          coords: { x: x2, y: y2 },
          rgb: { r, g, b },
        } = e.data;
        image.data[4 * (y2 * 100 + x2)] = r;
        image.data[4 * (y2 * 100 + x2) + 1] = g;
        image.data[4 * (y2 * 100 + x2) + 2] = b;
        image.data[4 * (y2 * 100 + x2) + 3] = 0xff;
        x++;
        if (x >= canvas.width) {
          y++;
          if (y === canvas.width) {
            done = true;
          }
          x = 0;
        }
      }
      if (!done) worker.postMessage([x, y]);
    };
    workers.push(worker);
  }
}

run();
