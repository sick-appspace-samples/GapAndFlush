------------------------------------------------------------------------
-- General settings
------------------------------------------------------------------------
-- Search line along gap, profiles extracted orthogonal to this line.
local lineStart = Point.create(256, 150)
local lineStop = Point.create(256, 800)

-- Probe length (mm, in each direction), probe line extends
-- this much on either side of the search line.
local probeDist = 10

-- Probe samples per mm
local samplesPerMM = 20

-- Profile extraction, number of profiles and step length between each
local stepLength = 5
local nP = 110

-- circular edge fitting parameters
local circleTargetRadius = 1.0
local inlierMargin = 0.10*circleTargetRadius

-- Load data
local cont = Image.readIconData('resources/0') --("C:\\Images\\ProfileImages\\0")
local range = cont:getObject('Ranger3Range.Range')
range = Image.crop(range, 0, 0, 2560, 950)

------------------------------------------------------------------------
-- Derived settings and view decorations
------------------------------------------------------------------------
-- number of samples to fit circle to (corresponding to radius of fitted circle)
local nSamplesCircle = math.ceil(circleTargetRadius*samplesPerMM)

-- General scale of decorations, based on the target circle radius
local decorationScale = circleTargetRadius

local v = View.create("view1")
local sDeco = View.ShapeDecoration.create()
sDeco:setLineWidth(50*decorationScale)

local v2 = View.create("view2")
local gDeco = View.GraphDecoration.create()
gDeco:setAspectRatio("EQUAL")
gDeco:setYBounds(-3, 0.5)

local tDeco = View.TextDecoration.create()
tDeco:setSize(1*decorationScale)
tDeco:setPosition(200, 300)

local tDeco2 = View.TextDecoration.create()
tDeco2:setSize(0.7*decorationScale)
tDeco2:setColor(200, 100, 0, 100)

local sDecoG = View.ShapeDecoration.create()
sDecoG:setLineWidth(inlierMargin)
sDecoG:setLineColor(200, 100, 0, 100)

local sDecoF = View.ShapeDecoration.create()
sDecoF:setLineWidth(inlierMargin)
sDecoF:setLineColor(200,0,200)

local tDecoF = View.TextDecoration.create()
tDecoF:setSize(0.7*decorationScale)
tDecoF:setColor(200,0,200)

local gDecoF = View.GraphDecoration.create()
gDecoF:setGraphColor(200,0,200)


------------------------------------------------------------------------
-- Functions
------------------------------------------------------------------------
---@param profile Profile
---@param nS int
---@param dir int
---@param targetRadius float
---@param outlierMargin float
---@return Shape
---@return float
---@return float
local function fitCircle(profile, nS, dir, targetRadius, outlierMargin)
  local circle
  local cQuality
  local poly

  local validVec = profile:getValidFlag()
  local validFlagProf = Profile.createFromVector(validVec)
  local validInd = Profile.findEqual(validFlagProf, 1.0, 0.1)
  local firstValid
  if dir  == 1 then
    firstValid = validInd[1] + 1
  else
    firstValid = validInd[#validInd] + 1
  end

  local profCirc
  if dir == 1 then
    profCirc = profile:crop(firstValid-1, firstValid-1 + nS - 1)
  else
    profCirc = profile:crop(firstValid-1 - nS + 1, firstValid-1)
  end

  local sf = Image.ShapeFitter.create()
  sf:setIterations(1000)
  sf:setOutlierMargin(outlierMargin, "ABSOLUTE")
  circle, cQuality = sf:fitCircle(profCirc:toPoints(), targetRadius, targetRadius)

  -- Fit polynomial
  local pLength = Profile.getSize(profile)
  local profPoly
  if dir == 1 then
    profPoly = profile:crop(firstValid-1 + nS, pLength-1)
  else
    profPoly = profile:crop(0, firstValid-1 - nS)
  end

  local CF = Profile.CurveFitter.create()
  poly = CF:fitPolynomial(profPoly, 1)

  return circle, cQuality, poly
end

------------------------------------------------------------------------
local function extractProfiles()
  local profiles = {}
  local probeLines = {}

  -- calculate direction vectors
  local lineDir = Point.subtract(lineStop, lineStart)
  lineDir = Point.normalize(lineDir)
  local probeDir = Point.create(lineDir:getY(), -lineDir:getX())

  -- calculate profile range value offsets to normalize extracted profiles
  -- fit a line to one side of the gap
  local pSC = Point.add(lineStart,Point.multiplyConstant(lineDir, 1*stepLength))
  local probeStart = Point.subtract(pSC, Point.multiplyConstant(probeDir, -probeDist))
  pSC = Point.add(lineStart,Point.multiplyConstant(lineDir, nP*stepLength))
  local probeStop = Point.subtract(pSC, Point.multiplyConstant(probeDir, -probeDist))
  local probeLine = Shape.createLineSegment(probeStart, probeStop)
  local profileAlong = Image.extractProfile(range, probeLine, nP)

  profileAlong:convertCoordinateType("IMPLICIT_1D")
  local CF = Profile.CurveFitter.create()
  CF:setFitMode("RANSAC")
  CF:setOutlierMargin(2, "ABSOLUTE")
  local lineA = CF:fitLine(profileAlong)
  local offsetA, slopeA = lineA:getLineParameters()

  -- Extract profile across the edge
  for k = 1, nP do
    pSC = Point.add(lineStart,Point.multiplyConstant(lineDir, k*stepLength))
    probeStart = Point.subtract(pSC, Point.multiplyConstant(probeDir, -probeDist))
    probeStop = Point.subtract(pSC, Point.multiplyConstant(probeDir, probeDist))
    probeLine = Shape.createLineSegment(probeStart, probeStop)

    local pTemp = Image.extractProfile(range, probeLine, 2*probeDist*samplesPerMM)
    pTemp:addConstantInplace(-(offsetA + k*stepLength*slopeA))
    profiles[#profiles + 1] = pTemp
    probeLines[#probeLines + 1] = probeLine
  end

  return profiles, probeLines
end

------------------------------------------------------------------------
-- Main
------------------------------------------------------------------------
Image.setMissingDataFlag(range, 1)
local profiles, probeLines = extractProfiles()

-- Collect neighboring profiles and calculate gap and flush on aggregated profiles.
local gapValues = {}
local aggDist = 4
for k = (1+aggDist), (#profiles - aggDist) do
  local profAggs = {}
  for kk = -aggDist, aggDist do
    profAggs[#profAggs + 1] = profiles[k + kk]
  end
  local profAgg, validProf = Profile.aggregate(profAggs, "MEAN")
  local invalidIndex = Profile.threshold(validProf, -1, 3)
  Profile.setValidFlag(profAgg, invalidIndex, 0)
  invalidIndex = Profile.threshold(profAgg, nil, -4)
  Profile.setValidFlag(profAgg, invalidIndex, 0)

  profAgg:convertCoordinateType("IMPLICIT_1D")

  -- Fit circle
  print(profAgg:toString())
  -- There are some divergent values near the gap that disturbs the fitting of the circle
  profAgg = Profile.removeDivergentValues(profAgg, 1, 5)
  local profRight = profAgg:crop(probeDist*samplesPerMM, 2*probeDist*samplesPerMM-1)
  local circleR, _, polyR = fitCircle(profRight, nSamplesCircle, 1, circleTargetRadius, inlierMargin)

  local profLeft = profAgg:crop(0, probeDist*samplesPerMM)
  local circleL, _, polyL = fitCircle(profLeft, nSamplesCircle, -1, circleTargetRadius, inlierMargin)

  -- Calculate gap
  local gap = -1
  local leftX
  local rightX
  if circleL ~= nil and circleR ~= nil then
    local pLeftC, radiusL = Shape.getCircleParameters(circleL)
    local pRightC, radiusR = Shape.getCircleParameters(circleR)

    leftX = pLeftC:getX() + radiusL
    rightX = pRightC:getX() - radiusR
    gap = rightX - leftX
  end
  local midwayPoint = (leftX+rightX)/2

  gapValues[k] = gap

  -- Calculate flush midway between the edges
  local polyParamR = Profile.Curve.getPolynomialParameters(polyR)
  local cHeightR = polyParamR[1] + polyParamR[2]*midwayPoint

  local polyParamL = Profile.Curve.getPolynomialParameters(polyL)
  local cHeightL = polyParamL[1] + polyParamL[2]*midwayPoint
  local flush = cHeightL - cHeightR

  v:clear()
  local id = v:addImage(range)
  v:addShape(lineStop,sDeco,nil,id)
  v:addShape(probeLines[k],sDeco,nil,id)
  v:addText(tostring(k),tDeco)
  v:present()

  v2:clear()
  local pId = v2:addProfile(profAgg, gDeco)
  v2:addShape(circleL, sDecoG)
  v2:addShape(circleR, sDecoG)
  if gap > 0 then
    -- Plot gap
    local ll = Shape.createLineSegment(Point.create(leftX,-3), Point.create(leftX,0))
    local rl = Shape.createLineSegment(Point.create(rightX,-3), Point.create(rightX,0))
    v2:addShape(ll, sDecoG)
    v2:addShape(rl, sDecoG)
    tDeco2:setPosition(leftX + 1, -2)
    local gs = string.format("%.2f", gap)
    local ts = "Gap: " .. gs .. " mm"
    v2:addText(ts, tDeco2)

    -- Plot flush
    local fL = Shape.createLineSegment(Point.create(midwayPoint, cHeightR),
                                       Point.create(midwayPoint, cHeightL))
    v2:addShape(fL, sDecoF)
    tDecoF:setPosition(leftX+1, 0.5)--(cHeightR + cHeightL)/2)
    local fs = string.format("%.2f", flush)
    local tsf = "Flush: " .. fs .. " mm"
    v2:addText(tsf, tDecoF)

    local profLeftRef = profAgg:crop(0, midwayPoint*samplesPerMM)
    local polyProfL = Profile.Curve.toProfile(polyL, profLeftRef)
    v2:addProfile(polyProfL,gDecoF,nil,pId)
    local profRightRef = profAgg:crop(midwayPoint*samplesPerMM, 2*probeDist*samplesPerMM-1)
    local polyProfR = Profile.Curve.toProfile(polyR, profRightRef)
    v2:addProfile(polyProfR,gDecoF,nil,pId)
  end
  v2:present()
end


-- Plot gap profile
local gapProfile = Profile.createFromVector(gapValues)
v:clear()
v:addProfile(gapProfile)
tDecoF:setPosition(33, 6.35)
tDecoF:setSize(5)
v:addText("Gap in mm", tDecoF)
v:present()


print("App finished")