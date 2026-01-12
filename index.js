const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");
let image = ctx.createImageData(100, 100);
let done;

function draw() {
  ctx.putImageData(image, 0, 0);
  if (!done) requestAnimationFrame(draw);
}

async function run() {
  done = false;
  requestAnimationFrame(draw);

  let nextX = 0;
  let nextY = 0;
  let workers = [];
  for (let i = 0; i < navigator.hardwareConcurrency; ++i) {
    const worker = new Worker("worker.js", { type: "module" });
    worker.onmessage = (e) => {
      if (e.data) {
        const {
          coords: { x, y },
          rgb: { r, g, b },
        } = e.data;
        image.data[4 * (y * 100 + x)] = r;
        image.data[4 * (y * 100 + x) + 1] = g;
        image.data[4 * (y * 100 + x) + 2] = b;
        image.data[4 * (y * 100 + x) + 3] = 0xff;
      }
      if (!done) {
        worker.postMessage([nextX, nextY]);
        nextX++;
        if (nextX >= canvas.width) {
          nextY++;
          if (nextY === canvas.width) {
            done = true;
          }
          nextX = 0;
        }
      }
    };
    workers.push(worker);
  }
}

run();
