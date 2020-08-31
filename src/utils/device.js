export function onMobileDevice() {
  // Explicitly allow iPads/tablets
  if (/Tablet|iPad/i.test(navigator.userAgent)) {
    return false;
  }

  return /Mobile|Android|webOS|iPhone|iPod|BlackBerry|BB|PlayBook|IEMobile|Windows Phone|Kindle|Silk|Opera Mini/i.test(
    navigator.userAgent
  );
}
