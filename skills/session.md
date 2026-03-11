---
name: session
description: Append a session entry to the current session_*.md file. Use when the user asks to log, record, or update the session, or at natural milestones during a session.
---

Append a summary of what happened in this session to the current session file.

## How to append

Run:
```
python scripts/session.py append "<text>"
```

The script handles rotation automatically — if the file is too old a new one is created.

## What to write

Append a session continuation block in this format:

```markdown
---

## Session Continuation (YYYY-MM-DD ~HH:xx)

### What Was Worked On
- One bullet per distinct task or change made
- Be specific: name the file, function, or component changed
- Include what was added, fixed, or removed

### Decisions Made
- Architectural or design choices agreed with the user
- Why one approach was chosen over another (brief rationale)
- Any rules or constraints established for future work

### Pain Points / Errors Found
- Bugs discovered (even if not fixed this session)
- Flaky tests, wrong expected values, broken scripts
- Write "None." if nothing found

### Next Steps Agreed On
- Numbered list of concrete next actions
- Each item should be actionable, not vague
- Include any items that need user approval before proceeding
```

## When to write

Append at natural milestones during a session. Do not wait until the end. Triggers include:
- A feature or bug fix is completed
- A key decision is made
- A significant discussion concludes that changes understanding of the system (e.g. architecture, tooling, workflow decisions) — even if no code was written

## Rules

- Only include things that actually happened — no speculation or aspirational items
- "Next Steps" must be things explicitly agreed on with the user, not inferred
- Keep bullets concise but specific enough to be useful cold-start context next session
- Do not duplicate content already written earlier in the same session file
