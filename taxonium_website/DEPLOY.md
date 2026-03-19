# Adding Trees and Deploying the Website

## Overview

The website is a Vite/React app. The list of available trees is defined in two places that must both be updated:

1. **`src/trees.json`** — metadata for each tree (URL, title, description)
2. **`src/App.jsx`** — the `SHOWCASE_PATHS` array, which controls which trees appear on the homepage

After editing either file, you must rebuild and redeploy.

---

## Adding a new tree

### 1. Add an entry to `src/trees.json`

Add a new key/value pair. The key is the URL path (e.g. `"flu/HA-H1"`). The value is:

```json
"flu/HA-H1": {
  "protoUrl": "https://taxonium.badelita.com/tree_HA_H1.jsonl.gz",
  "title": "Flu HA H1",
  "description": "Influenza A HA subtype H1 tree.",
  "icon": "/assets/usher.png",
  "maintainerMessage": "Maintained by Vlad Badelita"
}
```

- `protoUrl` — direct URL to the `.jsonl.gz` file on your server
- `title` — shown in the tree list
- `description` — shown beneath the title
- `icon` — leave as `/assets/usher.png`
- `maintainerMessage` — short credit line

### 2. Add the key to `SHOWCASE_PATHS` in `src/App.jsx`

Find the `SHOWCASE_PATHS` array near the top of the file and add the new key in the appropriate position:

```js
const SHOWCASE_PATHS = [
  "flu/HA-H1",
  "flu/HA-H1-allseg",   // <-- add new entries here
  ...
];
```

Trees not listed in `SHOWCASE_PATHS` are still accessible by direct URL but won't appear on the homepage.

---

## Rebuild and deploy

Use the script below, or run the steps manually.

### Manual steps

```bash
# From the repo root
cd taxonium_website
npm run build
rsync -avzP dist/ vlad@raven:~/taxonium_website_dist/
```

Then on raven:

```bash
sudo cp -r ~/taxonium_website_dist/* /var/www/taxonium/
```

### Deploy script

Run `./deploy.sh` from the `taxonium_website/` directory (see `deploy.sh` in this folder).

---

## Checklist

- [ ] Tree file uploaded to `/var/www/taxonium/` on raven
- [ ] Entry added to `src/trees.json`
- [ ] Key added to `SHOWCASE_PATHS` in `src/App.jsx`
- [ ] `npm run build` run successfully
- [ ] `deploy.sh` run (or manual rsync + copy)
