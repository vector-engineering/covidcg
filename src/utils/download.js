// Download utilities

export function downloadBlobURL(blob_url, filename) {
  let link = window.document.createElement('a');
  link.setAttribute('href', blob_url);
  link.setAttribute('download', filename);
  link.style.visibility = 'hidden';
  window.document.body.appendChild(link);
  link.click();
  window.document.body.removeChild(link);
}
