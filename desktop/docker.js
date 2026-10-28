// Run the docker command line. Apps started from the desktop do not get the
// shell's PATH, so the usual Docker Desktop locations are added.

const { spawn } = require('child_process');
const fs = require('fs');
const path = require('path');

const EXTRA_PATHS = {
  win32: ['C:\\Program Files\\Docker\\Docker\\resources\\bin'],
  darwin: ['/usr/local/bin', '/opt/homebrew/bin', '/Applications/Docker.app/Contents/Resources/bin'],
  linux: ['/usr/bin', '/usr/local/bin', '/snap/bin'],
};

const env = {
  ...process.env,
  PATH: [...(EXTRA_PATHS[process.platform] || []), process.env.PATH].join(path.delimiter),
};

// Resolves to { code, out }; code -1 means docker could not be started at all.
function docker(args, onLine) {
  return new Promise((resolve) => {
    let out = '';
    const child = spawn('docker', args, { env, windowsHide: true });
    const take = (chunk) => {
      const text = chunk.toString();
      out += text;
      if (onLine) text.split(/\r?\n|\r/).filter(Boolean).forEach(onLine);
    };
    child.stdout.on('data', take);
    child.stderr.on('data', take);
    child.on('error', (error) => resolve({ code: -1, out: String(error) }));
    child.on('close', (code) => resolve({ code, out: out.trim() }));
  });
}

// Docker Desktop is installed but not running: start it. Returns false if there is nothing to start.
function startDockerDesktop() {
  if (process.platform === 'win32') {
    const exe = 'C:\\Program Files\\Docker\\Docker\\Docker Desktop.exe';
    if (!fs.existsSync(exe)) return false;
    spawn(exe, [], { detached: true, stdio: 'ignore' }).unref();
    return true;
  }
  if (process.platform === 'darwin') {
    if (!fs.existsSync('/Applications/Docker.app')) return false;
    spawn('open', ['-a', 'Docker'], { detached: true, stdio: 'ignore' }).unref();
    return true;
  }
  return false;
}

module.exports = { docker, startDockerDesktop };
