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

  const worker = new Worker("worker.js", { type: "module" });
  let x = 0;
  let y = 0;
  worker.onmessage = (e) => {
    if (e.data) {
      const { r, g, b } = e.data;
      image.data[4 * (y * 100 + x)] = Math.floor(r * 255);
      image.data[4 * (y * 100 + x) + 1] = Math.floor(g * 255);
      image.data[4 * (y * 100 + x) + 2] = Math.floor(b * 255);
      image.data[4 * (y * 100 + x) + 3] = 0xff;
      x++;
      if (x >= canvas.width) {
        y++;
        x = 0;
      }
    }
    if (y < canvas.height) worker.postMessage([x, y]);
    else done = true;
  };
}

run();
