// ─────────────────────────────────────────────────────────────────────
//  Julia for HPC site — shared script
//  Drives: theme toggle, mobile sidebar, copy buttons, in-page
//          scrollspy, and cross-page active-link highlight.
// ─────────────────────────────────────────────────────────────────────

// ── Build info (commit hash + last-updated stamp in the footer) ────
// Populates any `<div id="build-info">` on the page from version.json.
// Silently no-ops on pages that don't have the element.
(function () {
  const el = document.getElementById('build-info');
  if (!el) return;
  fetch('version.json', { cache: 'no-store' })
    .then(r => r.ok ? r.json() : Promise.reject(r.status))
    .then(v => {
      const started = new Date(v.started_at).toLocaleString();
      el.textContent = `commit ${v.commit_short} · last updated ${started}`;
    })
    .catch(() => { el.textContent = 'build info unavailable'; });
})();

// ── Theme toggle ────────────────────────────────────────────────────
// The no-flash initial application lives in an inline <script> in each
// page's <head> (it must run before first paint). This handler just
// wires the click that flips the data-theme attribute and persists.
(function () {
  const btn = document.getElementById('theme-toggle');
  if (!btn) return;
  btn.addEventListener('click', () => {
    const isDark = document.documentElement.getAttribute('data-theme') === 'dark';
    if (isDark) {
      document.documentElement.removeAttribute('data-theme');
      localStorage.setItem('theme', 'light');
    } else {
      document.documentElement.setAttribute('data-theme', 'dark');
      localStorage.setItem('theme', 'dark');
    }
  });
})();

// ── Mobile sidebar toggle ────────────────────────────────────────────
(function () {
  const toggle  = document.getElementById('menu-toggle');
  const sidebar = document.getElementById('sidebar');
  const overlay = document.getElementById('sidebar-overlay');
  if (!toggle || !sidebar) return;

  const openMenu  = () => { sidebar.classList.add('open');    toggle.setAttribute('aria-expanded', 'true');  if (overlay) overlay.style.display = 'block'; };
  const closeMenu = () => { sidebar.classList.remove('open'); toggle.setAttribute('aria-expanded', 'false'); if (overlay) overlay.style.display = ''; };

  toggle.addEventListener('click', () => sidebar.classList.contains('open') ? closeMenu() : openMenu());
  if (overlay) overlay.addEventListener('click', closeMenu);

  sidebar.querySelectorAll('nav a').forEach(a => {
    a.addEventListener('click', () => { if (window.innerWidth <= 820) closeMenu(); });
  });
})();

// ── Copy-to-clipboard ────────────────────────────────────────────────
document.querySelectorAll('.copy-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    const block = btn.closest('.code-block');
    const panel = btn.closest('.deriv-code-panel');
    const text  = block ? block.querySelector('pre code')?.innerText
                : panel ? panel.querySelector('pre')?.innerText
                : null;
    if (!text) return;
    navigator.clipboard.writeText(text).then(() => {
      const label = btn.textContent;
      btn.textContent = 'Copied!';
      btn.classList.add('copied');
      setTimeout(() => { btn.textContent = label || 'Copy'; btn.classList.remove('copied'); }, 1800);
    });
  });
});

// ── Cross-page active link ───────────────────────────────────────────
// Marks the sidebar link whose href matches this page's filename.
// Falls back to comparing the basename so /pde-heat.html and pde-heat.html
// both resolve. Runs before scrollspy so the in-page tracker can refine it
// when there are multiple #anchors on the same page.
(function () {
  const sidebar = document.getElementById('sidebar');
  if (!sidebar) return;
  const here = (location.pathname.split('/').pop() || 'index.html').toLowerCase();

  sidebar.querySelectorAll('nav a[href]').forEach(a => {
    const href = a.getAttribute('href');
    if (!href || href.startsWith('#')) return;          // in-page anchors handled below
    // Strip query/hash to compare basenames cleanly
    const target = href.split('#')[0].split('?')[0].split('/').pop().toLowerCase();
    if (target === here) a.classList.add('active');
  });
})();

// ── In-page scrollspy (for chapter pages that have many h2/h3 anchors) ──
// Only enabled when the sidebar contains href="#…" links that match
// ids in the document. Otherwise no-op.
(function () {
  const sidebar      = document.getElementById('sidebar');
  if (!sidebar) return;
  const anchorLinks  = Array.from(sidebar.querySelectorAll('nav a[href^="#"]'));
  if (anchorLinks.length === 0) return;

  const linkFor = new Map();
  const targets = [];
  for (const a of anchorLinks) {
    const id = a.getAttribute('href').slice(1);
    const el = document.getElementById(id);
    if (el) { linkFor.set(id, a); targets.push(el); }
  }
  if (targets.length === 0) return;

  let lastActive = null;
  function update() {
    const triggerY = Math.max(80, window.innerHeight * 0.25);
    let current = targets[0];
    for (const el of targets) {
      if (el.getBoundingClientRect().top <= triggerY) current = el;
      else break;
    }
    const link = linkFor.get(current.id);
    if (link === lastActive) return;
    anchorLinks.forEach(a => a.classList.remove('active'));
    if (link) {
      link.classList.add('active');
      lastActive = link;
      // Keep the active link visible inside the sidebar
      const lr = link.getBoundingClientRect();
      const sr = sidebar.getBoundingClientRect();
      const margin = 40;
      if (lr.top < sr.top + margin)       sidebar.scrollTop += lr.top - sr.top - margin;
      else if (lr.bottom > sr.bottom - margin) sidebar.scrollTop += lr.bottom - sr.bottom + margin;
    }
  }
  let ticking = false;
  const onScroll = () => {
    if (ticking) return;
    ticking = true;
    requestAnimationFrame(() => { update(); ticking = false; });
  };
  window.addEventListener('scroll', onScroll, { passive: true });
  window.addEventListener('resize', onScroll, { passive: true });
  window.addEventListener('hashchange', update);
  window.addEventListener('load', update);
  update();
})();
