// AromaCards hors connexion : réseau d'abord (la page reste à jour), copie locale sans réseau ou si le réseau traîne.
// Le nom du cache porte la version publiée : à chaque publication, le nouveau service worker remplace l'ancien
// et supprime son cache (personne ne reste bloqué sur une ancienne copie).
const CACHE = "aromacards-05e4030db72c";
const PRECACHE = ["./", "index.html", "manifest.webmanifest", "favicon.svg", "favicon-32.png", "apple-touch-icon.png",
                  "icon-192.png", "icon-512.png", "icon-maskable-192.png", "icon-maskable-512.png"];
self.addEventListener("install", e => {
  self.skipWaiting();
  e.waitUntil(caches.open(CACHE).then(c => c.addAll(PRECACHE.map(u => new Request(u, {cache: "reload"})))));
});
self.addEventListener("activate", e => e.waitUntil(
  caches.keys().then(keys => Promise.all(keys.filter(k => k !== CACHE).map(k => caches.delete(k)))).then(() => self.clients.claim())));
self.addEventListener("fetch", e => {
  const url = new URL(e.request.url);
  if (e.request.method !== "GET" || url.origin !== location.origin) return;
  if (url.pathname.endsWith("/version.json")) return;  // contrôle de version : toujours le réseau, jamais de copie
  const nav = e.request.mode === "navigate";
  // pages : revalidées auprès du serveur (pas de copie HTTP périmée)
  const net = fetch(nav ? new Request(url.origin + url.pathname, {cache: "no-cache", credentials: "same-origin"}) : e.request)
    .then(r => { if (r.ok){ const copy = r.clone(); caches.open(CACHE).then(c => c.put(nav ? url.origin + url.pathname : e.request, copy)); } return r; });
  const local = () => caches.match(e.request, {ignoreSearch: true}).then(r => r || (nav ? caches.match("index.html") : undefined));
  if (!nav){ e.respondWith(net.catch(() => local().then(r => r || Response.error()))); return; }
  // page : au-delà de 4 s sans réponse, la copie locale s'affiche (un bandeau proposera la mise à jour)
  e.respondWith(new Promise(resolve => {
    let done = false;
    const give = r => { if (r && !done){ done = true; resolve(r); } };
    const timer = setTimeout(() => local().then(give), 4000);
    net.then(r => { clearTimeout(timer); give(r); })
       .catch(() => local().then(r => { clearTimeout(timer); give(r || Response.error()); }));
  }));
});
