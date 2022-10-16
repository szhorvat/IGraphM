
-- Generates separate table of contents for each section and susbection.
-- Adds links from sub-headings to the parent heading.
-- Based on https://groups.google.com/d/msg/pandoc-discuss/YA20DfSuQdY/fLZkC4P7AQAJ
-- Assumes that each heading has a unique ID.
-- Requires Pandoc 2.8 or later to ensure that headings are traversed in order.

local headings = {{}, {}, {}, {}, {}} -- initialize to store multiple levels
local parents = {}
local current_chapter = nil

-- Collects headings of 'chapter_level' into 'headings'
local function collect_headings (chapter_level)
  return function (head)
    if head.level <= chapter_level then
      local id = head.identifier
      current_chapter = {
        chapter = id,
        toc = {},
      }
      headings[chapter_level][id] = current_chapter
    elseif head.level == chapter_level + 1 then
      if current_chapter then
        local toc = current_chapter.toc
        toc[#toc+1] = head
        parents[head.identifier] = current_chapter.chapter
      end
    end
    return nil
  end
end

local function build_toc (heads, level)
  local toc = {}
  for _, head in ipairs(heads) do
    local entry = {
      pandoc.Plain{
        pandoc.Link(
          head.content:clone(), -- text
          '#' .. head.identifier, -- target
          "", -- empty title
          pandoc.Attr(
            "", -- empty identifier
            {'local-toc-link-' .. level} -- class
          )
        )
      }
    }
    toc[#toc+1] = entry
  end
  return pandoc.Div(
    { pandoc.BulletList(toc) },
    pandoc.Attr( "", {'local-toc-' .. level} )
  )
end

-- Insert table of contents for headings at 'chapter_level'
-- Must be run after 'collect_headings'
local function insert_toc (chapter_level)
  return function (head)
    if head.level <= chapter_level then
      local id = head.identifier
      if headings[chapter_level][id] then
        local toc = build_toc(
          headings[chapter_level][id].toc,
          chapter_level
        )
        return {head, toc}
      end
    end
    return nil
  end
end


-- Add a link back to the parent section for each section heading
local function head_add_backlink (head) 
  if parents[head.identifier] then
    head.content:extend(
      {
        pandoc.Link(
          pandoc.Span("▴", {class = 'link-back-triangle'}),
          "#" .. parents[head.identifier],
          "",
          {class = 'link-back'}
        )
      }
    )
  end
  return head
end


return {  
  { Header = collect_headings(1) },
  { Header = insert_toc(1) },
  { Header = collect_headings(2) },
  { Header = insert_toc(2) },
  { Header = head_add_backlink },
}
