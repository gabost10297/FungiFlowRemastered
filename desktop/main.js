// FungiFlow desktop app. It does not run the pipeline itself: it prepares a
// FungiFlow folder in Documents, makes sure the Docker images are there, starts
// the Shiny app container (image fungiflow-viz) and shows it in a window.
//
// Containers:
//   fungiflow-data   created, never started; holds the folder (/data) and Docker
//                    socket mounts that every other container takes with --volumes-from
//   fungiflow-app    the Shiny page; removed when the window closes
//   fungiflow_runner_<run>, fungiflow_<run>_<step>
//                    a pipeline run, started from the page; keeps going after the window closes

const { app, BrowserWindow, Menu, shell, dialog, ipcMain } = require('electron');
const fs = require('fs');
const http = require('http');
const net = require('net');
const path = require('path');
const { docker, startDockerDesktop } = require('./docker');
const settings = require('./package.json').fungiflow;

const VERSION = app.getVersion();
const IMAGES = ['fungiflow-remastered', 'fungiflow-ngspeciesid', 'fungiflow-viz'];
const APP_CONTAINER = 'fungiflow-app';
const DATA_CONTAINER = 'fungiflow-data';
const DOCKER_DOWNLOAD = 'https://www.docker.com/products/docker-desktop/';
const REPOSITORY = 'https://github.com/gabost10297/FungiFlowRemastered';
// FUNGIFLOW_LOCAL_IMAGES=1 uses images built with "docker build" instead of downloading them.
const LOCAL_IMAGES = process.env.FUNGIFLOW_LOCAL_IMAGES === '1';

// Documents, unless OneDrive syncs it: the databases and results are large and the
// run folders contain symlinks, so then the home folder is used.
const documents = app.getPath('documents');
const dataDir = path.join(/onedrive/i.test(documents) ? app.getPath('home') : documents, 'FungiFlow');
let win = null;

function status(text, kind = 'info') {
  if (kind !== 'detail') console.log(`[${kind}] ${text}`);
  if (win && !win.isDestroyed()) win.webContents.send('status', { text, kind });
}

class Problem extends Error {
  constructor(message, kind = 'error') {
    super(message);
    this.kind = kind;
  }
}

const sleep = (ms) => new Promise((done) => setTimeout(done, ms));

async function ensureDocker() {
  const ready = async () => (await docker(['info', '--format', '{{.ServerVersion}}'])).code === 0;
  if (await ready()) return;

  const found = await docker(['--version']);
  if (found.code === -1) {
    throw new Problem('Docker Desktop is not installed. FungiFlow needs it to run the ' +
      'pipeline tools. Install it, start it once, then press Try again.', 'docker');
  }
  if (startDockerDesktop()) status('Starting Docker Desktop. This can take a minute...');
  else status('Waiting for the Docker engine...');
  for (let i = 0; i < 60; i += 1) {
    await sleep(3000);
    if (await ready()) return;
  }
  throw new Problem('Docker is installed but not running. Start Docker Desktop ' +
    '(on Linux: the docker service), wait until it says it is running, then press Try again.', 'docker');
}

// Copies the scripts shipped with the app into the FungiFlow folder.
// config.sh and the UNITE archive are copied only once, so edits and downloads stay.
function projectSource() {
  return app.isPackaged ? path.join(process.resourcesPath, 'project') : path.join(__dirname, '..');
}

function copyTree(from, to, keep = () => true) {
  const info = fs.statSync(from);
  if (info.isDirectory()) {
    fs.mkdirSync(to, { recursive: true });
    for (const name of fs.readdirSync(from)) copyTree(path.join(from, name), path.join(to, name), keep);
    return;
  }
  if (!keep(from)) return;
  if (/\.(sh|py|R)$/.test(from)) {
    // Scripts run in Linux containers; a Windows checkout may have added \r.
    fs.writeFileSync(to, fs.readFileSync(from, 'utf8').replace(/\r\n/g, '\n'));
  } else {
    fs.copyFileSync(from, to);
  }
}

function ensureProject() {
  const source = projectSource();
  const stamp = path.join(dataDir, '.fungiflow-version');
  fs.mkdirSync(path.join(dataDir, 'runs'), { recursive: true });

  const current = fs.existsSync(stamp) && fs.readFileSync(stamp, 'utf8').trim() === VERSION;
  if (!current || !app.isPackaged) {
    fs.rmSync(path.join(dataDir, 'scripts'), { recursive: true, force: true });
    copyTree(path.join(source, 'scripts'), path.join(dataDir, 'scripts'),
      (file) => !file.includes('__pycache__'));
    copyTree(path.join(source, 'viz'), path.join(dataDir, 'viz'), (file) => file.endsWith('.R'));
    copyTree(path.join(source, 'run_pipeline.sh'), path.join(dataDir, 'run_pipeline.sh'));
    copyTree(path.join(source, 'README.md'), path.join(dataDir, 'README.md'));
  }
  if (!fs.existsSync(path.join(dataDir, 'config.sh'))) {
    copyTree(path.join(source, 'config.sh'), path.join(dataDir, 'config.sh'));
  }
  const unite = path.join(dataDir, 'database', 'unite');
  fs.mkdirSync(unite, { recursive: true });
  if (!fs.readdirSync(unite).some((name) => name.endsWith('.tgz'))) {
    copyTree(path.join(source, 'database', 'unite'), unite, (file) => file.endsWith('.tgz'));
  }
  fs.writeFileSync(stamp, `${VERSION}\n`);
}

// Images are published as ghcr.io/<owner>/<name>:<app version> and tagged with
// the short names that run_pipeline.sh and the page use.
async function ensureImages() {
  for (const [i, name] of IMAGES.entries()) {
    if (LOCAL_IMAGES) {
      if ((await docker(['image', 'inspect', name])).code !== 0) {
        throw new Problem(`FUNGIFLOW_LOCAL_IMAGES is set but there is no local image ${name}.`);
      }
      continue;
    }
    const remote = `${settings.registry}/${name}:${VERSION}`;
    if ((await docker(['image', 'inspect', remote])).code !== 0) {
      status(`Downloading the tools (${i + 1} of ${IMAGES.length}: ${name}). ` +
        'This happens once and can take a while (about 9 GB together).');
      // Only amd64 images exist; Apple Silicon Macs run them through emulation.
      const pull = await docker(['pull', '--platform', 'linux/amd64', remote], (line) => status(line, 'detail'));
      if (pull.code !== 0) {
        throw new Problem(`Could not download ${remote}. Check the internet connection, then press Try again.\n` +
          pull.out.split('\n').slice(-2).join(' '));
      }
    }
    await docker(['tag', remote, name]);
  }
}

// On Windows and macOS Docker Desktop writes files as the logged-in user whatever
// the container user is, and its socket belongs to root, so the containers run as root.
// On Linux they run as the user, plus the socket's group.
async function identity() {
  if (process.platform !== 'linux') return { user: '0:0', group: '' };
  const user = `${process.getuid()}:${process.getgid()}`;
  const socket = await docker(['run', '--rm', '-v', '/var/run/docker.sock:/var/run/docker.sock',
    'fungiflow-viz', 'stat', '-c', '%g', '/var/run/docker.sock']);
  return { user, group: socket.code === 0 ? socket.out.trim() : '' };
}

async function runningRuns() {
  const result = await docker(['ps', '--format', '{{.Names}}']);
  return result.code === 0 ? result.out.split('\n').filter((name) => name.startsWith('fungiflow_runner_')) : [];
}

const samePath = (a, b) => {
  const norm = (p) => p.replace(/\\/g, '/').replace(/\/+$/, '');
  return process.platform === 'linux' ? norm(a) === norm(b) : norm(a).toLowerCase() === norm(b).toLowerCase();
};

async function ensureDataContainer() {
  const mounted = await docker(['inspect', DATA_CONTAINER, '--format',
    '{{range .Mounts}}{{if eq .Destination "/data"}}{{.Source}}{{end}}{{end}}']);
  if (mounted.code === 0 && samePath(mounted.out, dataDir)) return;
  if (mounted.code === 0) {
    if ((await runningRuns()).length) {
      throw new Problem('A run from another FungiFlow folder is still going. Wait for it to finish, then press Try again.');
    }
    await docker(['rm', '-f', DATA_CONTAINER]);
  }
  const created = await docker(['create', '--name', DATA_CONTAINER,
    '-v', `${dataDir}:/data`, '-v', '/var/run/docker.sock:/var/run/docker.sock',
    'fungiflow-viz', 'true']);
  if (created.code !== 0) throw new Problem(`Docker could not use ${dataDir}:\n${created.out}`);
}

function freePort(first) {
  return new Promise((resolve, reject) => {
    const attempt = (port) => {
      if (port > first + 50) return reject(new Problem('No free network port for the page.'));
      const server = net.createServer();
      server.once('error', () => attempt(port + 1));
      server.once('listening', () => server.close(() => resolve(port)));
      server.listen(port, '127.0.0.1');
    };
    attempt(first);
  });
}

function answers(url) {
  return new Promise((resolve) => {
    http.get(url, (response) => { response.resume(); resolve(response.statusCode === 200); })
      .on('error', () => resolve(false));
  });
}

async function startApp(id) {
  await docker(['rm', '-f', APP_CONTAINER]);
  const port = await freePort(settings.port);
  const args = ['run', '-d', '--name', APP_CONTAINER, '--volumes-from', DATA_CONTAINER,
    '-p', `127.0.0.1:${port}:3838`, '--user', id.user,
    '-e', 'HOME=/tmp', '-e', `RUN_AS=${id.user}`, '-e', `DATA_CONTAINER=${DATA_CONTAINER}`];
  if (id.group) args.push('--group-add', id.group, '-e', `DOCKER_GID=${id.group}`);
  args.push('fungiflow-viz');
  const started = await docker(args);
  if (started.code !== 0) throw new Problem(`The page did not start:\n${started.out}`);

  const url = `http://127.0.0.1:${port}`;
  for (let i = 0; i < 90; i += 1) {
    if (await answers(url)) return url;
    await sleep(1000);
  }
  const logs = await docker(['logs', '--tail', '5', APP_CONTAINER]);
  throw new Problem(`The page did not answer:\n${logs.out}`);
}

async function launch() {
  try {
    status('Checking Docker...');
    await ensureDocker();
    status(`Preparing ${dataDir}...`);
    ensureProject();
    await ensureImages();
    const id = await identity();
    await ensureDataContainer();
    status('Starting FungiFlow...');
    const url = await startApp(id);
    await win.loadURL(url);
  } catch (error) {
    status(error.message, error.kind || 'error');
  }
}

function showStatusPage() {
  return win.loadFile(path.join(__dirname, 'status.html')).then(launch);
}

function buildMenu() {
  const template = [
    ...(process.platform === 'darwin' ? [{ role: 'appMenu' }] : []),
    {
      label: 'File',
      submenu: [
        { label: 'Open FungiFlow folder', click: () => shell.openPath(dataDir) },
        { label: 'Restart FungiFlow', click: () => showStatusPage() },
        { type: 'separator' },
        { role: 'quit' },
      ],
    },
    { role: 'editMenu' },
    {
      label: 'View',
      submenu: [{ role: 'reload' }, { role: 'toggleDevTools' }, { type: 'separator' },
        { role: 'resetZoom' }, { role: 'zoomIn' }, { role: 'zoomOut' }],
    },
    {
      label: 'Help',
      submenu: [
        { label: 'Read me', click: () => shell.openPath(path.join(dataDir, 'README.md')) },
        { label: 'FungiFlow on GitHub', click: () => shell.openExternal(REPOSITORY) },
        { label: `Version ${VERSION}`, enabled: false },
      ],
    },
  ];
  Menu.setApplicationMenu(Menu.buildFromTemplate(template));
}

function createWindow() {
  win = new BrowserWindow({
    width: 1400,
    height: 950,
    title: 'FungiFlow',
    webPreferences: { preload: path.join(__dirname, 'preload.js'), contextIsolation: true },
  });
  win.webContents.setWindowOpenHandler(({ url }) => {
    shell.openExternal(url);
    return { action: 'deny' };
  });
  showStatusPage();
}

if (!app.requestSingleInstanceLock()) {
  app.quit();
} else {
  app.on('second-instance', () => {
    if (win) {
      if (win.isMinimized()) win.restore();
      win.focus();
    }
  });

  ipcMain.on('retry', () => showStatusPage());
  ipcMain.on('download-docker', () => shell.openExternal(DOCKER_DOWNLOAD));
  ipcMain.on('open-folder', () => {
    fs.mkdirSync(dataDir, { recursive: true });
    shell.openPath(dataDir);
  });

  app.whenReady().then(() => {
    buildMenu();
    createWindow();
  });

  let stopped = false;
  app.on('before-quit', async (event) => {
    if (stopped) return;
    event.preventDefault();
    if ((await runningRuns()).length) {
      dialog.showMessageBoxSync({
        type: 'info',
        message: 'A pipeline run is still going.',
        detail: 'It continues in the background while Docker is running. Open FungiFlow again to follow it.',
      });
    }
    await docker(['rm', '-f', APP_CONTAINER]);
    stopped = true;
    app.quit();
  });

  app.on('window-all-closed', () => app.quit());
}
