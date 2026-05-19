// ── Mobile sidebar toggle ─────────────────────────────────────────────────
(function () {
  const toggle  = document.getElementById('menu-toggle');
  const sidebar = document.getElementById('sidebar');
  const overlay = document.getElementById('sidebar-overlay');
  if (!toggle || !sidebar) return;

  function openMenu() {
    sidebar.classList.add('open');
    toggle.setAttribute('aria-expanded', 'true');
    if (overlay) overlay.style.display = 'block';
  }

  function closeMenu() {
    sidebar.classList.remove('open');
    toggle.setAttribute('aria-expanded', 'false');
    if (overlay) overlay.style.display = '';
  }

  toggle.addEventListener('click', () => {
    sidebar.classList.contains('open') ? closeMenu() : openMenu();
  });

  if (overlay) overlay.addEventListener('click', closeMenu);

  // Close when a nav link is tapped on mobile
  sidebar.querySelectorAll('nav a').forEach(a => {
    a.addEventListener('click', () => {
      if (window.innerWidth <= 820) closeMenu();
    });
  });
})();

// Copy-to-clipboard — works for both .code-block and .deriv-code-panel
document.querySelectorAll('.copy-btn').forEach(btn => {
  btn.addEventListener('click', () => {
    let text;
    const block = btn.closest('.code-block');
    const panel = btn.closest('.deriv-code-panel');

    if (block) {
      text = block.querySelector('pre code').innerText;
    } else if (panel) {
      text = panel.querySelector('pre').innerText;
    }
    if (!text) return;

    navigator.clipboard.writeText(text).then(() => {
      btn.textContent = 'Copied!';
      btn.classList.add('copied');
      setTimeout(() => { btn.textContent = 'Copy'; btn.classList.remove('copied'); }, 2000);
    });
  });
});

// ── Highlight active sidebar link based on scroll position ─────────────────
// Only tracks elements that have a corresponding sidebar link; picks the
// topmost section whose top edge is above the trigger line.
(function () {
  const sidebar      = document.getElementById('sidebar');
  const sidebarLinks = Array.from(document.querySelectorAll('#sidebar nav a[href^="#"]'));
  if (!sidebar || sidebarLinks.length === 0) return;

  // Map: section id → its anchor link (only ids that exist in the DOM)
  const linkFor = new Map();
  const targets = [];
  for (const a of sidebarLinks) {
    const id = a.getAttribute('href').slice(1);
    const el = document.getElementById(id);
    if (el) {
      linkFor.set(id, a);
      targets.push(el);
    }
  }
  if (targets.length === 0) return;

  let lastActive = null;

  function update() {
    // Trigger line ~25% down the viewport — feels natural when reading
    const triggerY = Math.max(80, window.innerHeight * 0.25);

    // Walk targets in document order; the last one whose top is at or above
    // the trigger line is the current section.
    let current = targets[0];
    for (const el of targets) {
      if (el.getBoundingClientRect().top <= triggerY) {
        current = el;
      } else {
        break;
      }
    }

    const link = linkFor.get(current.id);
    if (link === lastActive) return;

    sidebarLinks.forEach(a => a.classList.remove('active'));
    if (link) {
      link.classList.add('active');
      lastActive = link;

      // Keep the active link visible inside the sidebar without disturbing
      // the main page scroll. Manual math beats scrollIntoView's variable
      // behaviour across browsers when the sidebar has overflow-y:auto.
      const lr = link.getBoundingClientRect();
      const sr = sidebar.getBoundingClientRect();
      const margin = 40;
      if (lr.top < sr.top + margin) {
        sidebar.scrollTop += lr.top - sr.top - margin;
      } else if (lr.bottom > sr.bottom - margin) {
        sidebar.scrollTop += lr.bottom - sr.bottom + margin;
      }
    }
  }

  let ticking = false;
  function onScroll() {
    if (ticking) return;
    ticking = true;
    requestAnimationFrame(() => { update(); ticking = false; });
  }

  window.addEventListener('scroll', onScroll, { passive: true });
  window.addEventListener('resize', onScroll, { passive: true });
  window.addEventListener('hashchange', update);
  window.addEventListener('load', update);
  update();
})();
