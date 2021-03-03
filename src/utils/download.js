// Download utilities

export function downloadBlobURL(blob_url, filename) {
  const link = window.document.getElementById('hidden-download-link');
  link.setAttribute('href', blob_url);
  link.setAttribute('download', filename);
  link.click();
}
