const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('fungiflow', {
  onStatus: (show) => ipcRenderer.on('status', (_event, status) => show(status)),
  retry: () => ipcRenderer.send('retry'),
  downloadDocker: () => ipcRenderer.send('download-docker'),
  openFolder: () => ipcRenderer.send('open-folder'),
});
