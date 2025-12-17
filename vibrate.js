const box = document.querySelector('#string-box')
const string = document.querySelector('#string')
// const svgWrapper = document.querySelector('.svg-wrapper');
const decay = 0.99;

const lineRect = box.getBoundingClientRect();
const svgX = 400 / (lineRect.right - lineRect.left);
const svgY = 200 / (lineRect.bottom - lineRect.top);

let lineY = 0;
let amplitude = 0;
let xPos = 0;
let isPlucking = false;
let time = 0.0;


document.addEventListener('mousemove', (e) => {
    const lineRect = box.getBoundingClientRect();
    lineY = lineRect.top + (lineRect.height / 2);

    mouseX = e.clientX;
    mouseY = e.clientY;

    const distanceY = mouseY - lineY;
    const absDistanceY = Math.abs(distanceY);
    const isWithinX = mouseX >= lineRect.left && mouseX <= lineRect.right;
    if (absDistanceY < 80 && isWithinX) {       // figure out if X is in the right area too 
        if (absDistanceY < 7) {
            isPlucking = true;
            console.log("a");
        }
        if (isPlucking) {
            xPos = mouseX - lineRect.left;
            amplitude = distanceY;
        }
    } else {
        isPlucking = false;
    }
})

function generateWavePath(curAmplitude) {
    const d = `M 0 100 `;

    const control = `Q ${xPos * svgX} ${100 + curAmplitude} `;

    const end = `400 100`;

    return d + control + end;
}

function animate() {

    if (isPlucking) {
        time = 0.0;
    } else {
        amplitude *= decay;
        time += 0.2;
    }

    const curAmplitude = Math.cos(time) * amplitude;

    string.setAttribute('d', generateWavePath(curAmplitude));

    requestAnimationFrame(animate); // need to cancel animationframe at some point
}

animate();
