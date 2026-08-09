# AGENTS.md

## Working agreements

### Baseline

1. Do not fabricate, manipulate any experimental data.
2. No permission to read and write any experimental data, not matter in what format. 
3. Result or knowledge provided must come with valid and last updated references.
4. Do not search things at the PC base directory. The only spaces you can search are the in the project folders and global `.codex` or `.claude` (such as skills and memeories)
5. Always as me for permission before you run, read, and edit any files in Nichol group box. 
6. You have no permission to edit things out side each project folder that the user are currently working on, unless the user permit you. 

### From NG lab manual

APS, NSF, DOD, and NIH all have guidelines about approperiate use of AI agent. The primary concerns are the following. 

1. An AI cannot be held responsible for any errors. As such, it cannot be redited with authorship. This means that any errors introduced through the use of AI, such as plagiarism or incorrect calculations or citations, are the responsibility of the human authors.
2. Any information you upload to an AI is effectively in the public domain. At a minimum, it means you cannot upload confidential documents like papers or proposals that you are reviewing, for example. However, it also means that the draft of the paper you just uploaded for edititing is no longer confidential.

The current guidance is the following:

- Substantial parts of the paper should not be uploaded to an AI agent before publication.
- If you use an AI for any part of the writing process, please let the PI know.
- Violations of these rules leading yo unintentional errors or plagiarism will be treated as academic ethics issue. 

#### The plagiarism rule

Plagiarism is never acceptable.

- Plagiarism means copying someone else’s work (words, pictures, thoughts, etc.) without providing
credit.
- When you cite a reference, that means you are quoting the ideas of the reference. It does not
mean that you are directly quoting words from the reference.
- If you must copy text, use quotation marks.
- Plagiarism applies to figures. 

## Coding agreement

1. Python (peps-skill) and Matlab coding (matlab-coding) should follow their own coding styling and guildline, see `skills\.coding`.
2. All the codes you write, should refer to the existing code as possible, especially for code for measurement, requiring knowledge inside you skill .ngskill.
3. All code you write you need to make sure the user understand completely the code. 
4. All code you write must be friendly to read. For exmaple, control the line spacing, no overstacking of lines and variables, etc. This belongs to the styling. 

## Conversation summaries and multi-agent task records

These rules apply to future host agents, subagents, test/verifier agents, and
single-agent conversations. Project instructions may add narrower requirements.
User authorization, confidentiality rules, experimental-data restrictions, and
filesystem permissions still take precedence.

### When to create a durable record

1. Do not create files for a small, self-contained conversation by default.
2. Create a durable record when the user requests one, work must survive context
   compaction or another session, multiple agents need shared context, the task
   has several dependent steps, or independent verification must be preserved.
3. A summary is a selective operational checkpoint, not a transcript. Include
   objectives, decisions, constraints, current state, authorized paths, task IDs,
   evidence, unresolved questions, and next steps. Exclude conversational noise.

### Location and directory layout

1. Store task records inside the current project, never in the global
   `~/.codex` or `~/.claude` directory. Global `.codex` or `.claude` is for Codex or Claude code instructions,
   configuration, skills, plugins, and reusable agent role definitions only.
2. Follow the project's existing documentation system of record. If
   `docs/exec-plans/` exists, use:

   ```text
   docs/exec-plans/
   ├── active/
   │   └── <YYYY-MM-DD-task-slug>/
   │       ├── plan.md
   │       ├── context/
   │       │   └── conversation-v001.md
   │       ├── tasks/
   │       │   └── T01-<task>.md
   │       ├── results/
   │       │   └── T01-<agent>-r01.md
   │       ├── verification/
   │       │   └── V01-<scope>.md
   │       └── final-summary.md
   └── completed/
   ```

3. If the project has no execution-plan convention, propose
   `docs/agent-runs/active/<YYYY-MM-DD-task-slug>/` with the same internal
   structure. Ask before creating it when existing permissions require approval.
4. Never place agent summaries or IO artifacts in experimental data, runtime
   metadata, generated binary, credential, cache, or human lab-notebook
   directories.
5. Reusable custom-agent definitions belong in project or global
   `.codex/agents/*.toml` or `.claude/agents/*.toml`, not in task-record Markdown.

### File ownership and single-writer rule

| Artifact | Sole editor | Purpose |
|---|---|---|
| `plan.md` | Host agent | Objective, task graph, assignments, status, decisions |
| `context/conversation-vNNN.md` | Host agent | Versioned shared context checkpoint |
| `tasks/TNN-*.md` | Host agent | Exact subtask, inputs, scope, and output contract |
| `results/TNN-*.md` | Assigned subagent | Independent subagent result and evidence |
| `verification/VNN-*.md` | Test/verifier agent | Independent checks, findings, and verdict |
| `final-summary.md` | Host agent | Final synthesis after verification |

1. One agent owns each file. Never have agents concurrently edit a shared file.
2. After a context, task, result, or verification file has been consumed, treat
   it as immutable. Create `conversation-v002.md`, `-r02.md`, or `V02-*.md`
   instead of silently rewriting history.
3. The host must not rewrite a subagent or verifier report. It may respond in
   `plan.md` or `final-summary.md`, or request a new revision.
4. If a subagent cannot write because it is read-only, it returns its structured
   result to the host. The host may persist it verbatim, only with authorization,
   and must identify the original author and the host as the persisting agent.

### Host-agent procedure

1. Confirm the user's objective, authorization, prohibited paths, and definition
   of done.
2. Create `plan.md` and the minimum versioned context needed for the run.
3. Divide work into bounded tasks with stable IDs (`T01`, `T02`, ...). Each task
   file names the context version, authorized inputs, prohibited paths, expected
   output file, acceptance criteria, and dependencies.
4. Give each subagent only the context and files needed for its task. Do not use
   delegation to widen permissions or bypass a user confirmation requirement.
5. Track assignments and state in `plan.md`; live messages are for coordination,
   while Markdown files are durable checkpoints.
6. Wait for required results, reconcile conflicts explicitly, and send the host's
   implementation or synthesis to an independent verifier when verification is
   required.
7. Write `final-summary.md` only after handling verifier findings. Preserve raw
   task and verification reports.
8. Move the complete run directory from `active/` to `completed/` only after the
   outcome, remaining limitations, and evidence state are recorded and the move
   is authorized.

### Subagent IO procedure

1. Read the assigned task file and its named context version first.
2. Stay inside the assigned scope. A host agent cannot grant authority the user
   did not grant, and a subagent inherits the parent session's restrictions.
3. Do not edit `plan.md`, shared context, another agent's task/result, or the
   final summary.
4. Write only the assigned `results/TNN-<agent>-rNN.md`, or return the same
   structure to the host when direct writes are unavailable.
5. A result records: task ID, context version, status, sources inspected, work
   performed, evidence, files changed (if authorized), checks run, limitations,
   unresolved questions, and recommended next action.
6. Send the host a short completion message with the result path and headline;
   keep the full evidence in the result artifact.

### Test/verifier-agent procedure

1. Be independent of the host's conclusion. Verify against the task acceptance
   criteria and authorized implementation/evidence, not only the host summary.
2. Remain read-only with respect to implementation files unless the user
   separately authorizes a fix. Write only the assigned verification report.
3. Record verification scope, environment/assumptions, commands or methods,
   observed evidence, pass/fail results, findings ordered by severity,
   reproduction details, untested areas, and a clear verdict.
4. Distinguish a confirmed defect from a hypothesis or missing evidence. Never
   claim experimental validation from static review or simulation.
5. Do not modify an earlier report after the host fixes an issue. Produce a new
   numbered verification report so the audit trail remains intact.

### Privacy, provenance, and permissions

1. Never put credentials, secrets, confidential text, restricted experimental
   data, measurement-derived fixtures, or copied full conversations into agent
   task records.
2. Summaries and results must retain source links, last-updated dates, evidence
   state, and human ownership where relevant.
3. Before reading or writing project task records, follow the project's approval
   rules. Direct-file IO is optional; use agent messaging when files are not
   authorized or persistence is unnecessary.
4. These procedures organize communication; they do not relax any safety,
   data-access, authorship, citation, or filesystem rule elsewhere in this file.


#### Others

1. When summarize, making plans, making prompts to subagents, IO with subagents, you need to write down markdown file in docs, following the rule for general agent development and AGENT.md. 
2. To save tokens, when a new agent load and run a lot shell code in terminal, that agent should summaize the results for the future agent to use that infomation to avoid repetitively calling information as waste tokens. 
3. Do not run all python file by agent. Python file should be run by human side and feedback will be send back to agent. So the agent need to give detialed instruction to the user that which python file need to run and what result they need to return. 